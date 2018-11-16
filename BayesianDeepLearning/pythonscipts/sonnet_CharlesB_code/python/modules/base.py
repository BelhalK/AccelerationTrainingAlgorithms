# Copyright 2017 The Sonnet Authors. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or  implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ============================================================================

"""Base class for TensorFlow snt.

This file contains the Abstract Base Class for defining Modules in TensorFlow.
A Module is an object that can be connected into the Graph multiple times
using the __call__ method, sharing variables automatically with no need to
explicitly use scopes or specify reuse=True.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import abc
import collections
import contextlib
import functools
import inspect
import weakref

# Dependency imports
import six
from sonnet.python.modules import base_info
from sonnet.python.modules import util
import tensorflow as tf

# Import error class from base_errors for backward compatibility.

from sonnet.python.modules.base_errors import Error
from sonnet.python.modules.base_errors import NotConnectedError
from sonnet.python.modules.base_errors import ParentNotBuiltError
from sonnet.python.modules.base_errors import IncompatibleShapeError
from sonnet.python.modules.base_errors import UnderspecifiedError
from sonnet.python.modules.base_errors import NotSupportedError
from sonnet.python.modules.base_errors import NotInitializedError
from sonnet.python.modules.base_errors import DifferentGraphError
from sonnet.python.modules.base_errors import ModuleInfoError
# pylint: enable=g-bad-import-order
# pylint: enable=unused-import


# Maps `tf.Graph` objects to a module call stack.
_MODULE_STACKS = weakref.WeakKeyDictionary()


# Maps `tf.Graph` objects to a stack of module connection observers.
_CONNECTION_OBSERVER_STACKS = weakref.WeakKeyDictionary()


@contextlib.contextmanager
def observe_connections(observer):
  """Notifies the observer whenever any Sonnet module is connected to the graph.

  If a module contains nested modules, the observer is notified once for each
  nested module, followed by the containing module.

  For example:

  ```python
  def logging_observer(connected_subgraph):
    logging.info(connected_subgraph.module.module_name)

  with snt.observe_connections(logging_observer):
    output = imagenet_module(input_tensor)
  ```

  Args:
    observer: Callable accepting a single argument. Will be called with a
    `ConnectedSubGraph` each time a module is connected to the graph.

  Yields:
    None: just yields control to the inner context.
  """
  connection_observer_stack = _CONNECTION_OBSERVER_STACKS.setdefault(
      tf.get_default_graph(), [])
  connection_observer_stack.append(observer)
  try:
    yield
  finally:
    connection_observer_stack.pop()



def _maybe_wrap_custom_getter(custom_getter, old_getter):
  """Wrap a call to a custom_getter to use the old_getter internally.

  Copied from [variable_scope._maybe_wrap_custom_getter](
  https://github.com/tensorflow/tensorflow/blob/master/tensorflow/python/
  ops/variable_scope.py#L1565)

  Args:
    custom_getter: The wrapping custom getter.
    old_getter: The wrapped custom getter.

  Returns:
    A new custom getter that calls `old_getter` and then `custom_getter`.
  """
  if old_getter is None:
    return custom_getter

  # The new custom_getter should call the old one
  def wrapped_custom_getter(getter, *args, **kwargs):
    # Call:
    #  custom_getter(
    #    lambda: old_getter(true_getter, ...), *args, **kwargs)
    # which means custom_getter will call old_getter, which
    # will call the true_getter, perform any intermediate
    # processing, and return the results to the current
    # getter, which will also perform additional processing.
    return custom_getter(functools.partial(old_getter, getter), *args, **kwargs)

  return wrapped_custom_getter


def _variable_tracking_custom_getter(getter, *args, **kwargs):
  """Custom getter that tracks variables created.

  This custom getter places any variables that `getter` creates into the
  `_all_variables` attribute of the `AbstractModule` that is on top of the
  module call stack. The module call stack is a graph-dependent stack that
  keeps track of the sonnet module call order.

  Note that this assumes that variables added appended to `tf.Graph`
  collections. This is a safe assumption to make because
  `tf.add_to_collection()` appends objects to collections, and `tf.Variable`
  uses `tf.add_to_collections()` to add itself to `tf.Graph` collections.

  Note that this assumes that all variables are added either the
  `tf.GraphKeys.GLOBAL_VARIABLES` or `tf.GraphKeys.LOCAL_VARIABLES` collection.

  Args:
    getter: The true getter or another custom getter.
    *args: See positional arguments for `tf.get_variable()`.
    **kwargs: See keyword arguments for `tf.get_variable()`.

  Returns:
    See docstring for `tf.get_variable()`.
  """
  # Get the module that is calling `tf.get_variable()`
  module_stack = _MODULE_STACKS[tf.get_default_graph()]
  module = module_stack[-1]

  # Get lists of local and global variables. We use `tf.get_collection_ref()`
  # instead of `tf.get_collection()` to avoid copying the collections.
  local_variables = tf.get_collection_ref(tf.GraphKeys.LOCAL_VARIABLES)
  global_variables = tf.get_collection_ref(tf.GraphKeys.GLOBAL_VARIABLES)

  num_local_vars_before = len(local_variables)
  num_global_vars_before = len(global_variables)

  out = getter(*args, **kwargs)

  # Add any local or global variables that have been created to `module`
  # pylint: disable=protected-access
  module._all_variables.update(local_variables[num_local_vars_before:])
  module._all_variables.update(global_variables[num_global_vars_before:])
  # pylint: enable=protected-access

  return out


@six.add_metaclass(abc.ABCMeta)
class AbstractModule(object):
  """Superclass for Sonnet Modules.

  This class defines the functionality that every module should implement,
  principally the `build` method which is wrapped using `tf.make_template`
  and called from `__call__`. Every time the module is called it will
  be connected into the graph but using the same shared set of variables, thanks
  to the template.

  For this to work correctly, the `build` implementation in the derived class
  must access all variables using `tf.get_variable`, not `tf.Variable`. The same
  set of variables must be created each time, if this is not the case an Error
  will be raised.

  Every subclass must call this class' `__init__` at the start of their
  `__init__`, passing the relevant name. If this step is omitted variable
  sharing will not work.
  """

  def __init__(self, _sentinel=None, custom_getter=None,
               name=None):  # pylint: disable=invalid-name
    """Performs the initialisation necessary for all AbstractModule instances.

    Every subclass of AbstractModule must begin their constructor with a call to
    this constructor, i.e.

    `super(MySubModule, self).__init__(custom_getter=custom_getter, name=name)`.

    If you instantiate sub-modules in __init__ you must create them within the
    `_enter_variable_scope` context manager to ensure they are in the module's
    variable scope. Alternatively, instantiate sub-modules in `_build`.

    Args:
      _sentinel: Variable that only carries a non-None value if `__init__` was
          called without named parameters. If this is the case, a deprecation
          warning is issued in form of a `ValueError`.
      custom_getter: Callable or dictionary of callables to use as
        custom getters inside the module. If a dictionary, the keys
        correspond to regexes to match variable names. See the `tf.get_variable`
        documentation for information about the custom_getter API.
      name: Name of this module. Used to construct the Templated build function.
          If `None` the module's class name is used (converted to snake case).

    Raises:
      TypeError: If `name` is not a string.
      TypeError: If a given `custom_getter` is not callable.
      ValueError: If `__init__` was called without named arguments.
    """
    if _sentinel is not None:
      raise ValueError("Calling AbstractModule.__init__ without named "
                       "arguments is not supported.")

    if name is None:
      name = util.to_snake_case(self.__class__.__name__)
    elif not isinstance(name, six.string_types):
      raise TypeError("Name must be a string, not {} of type {}.".format(
          name, type(name)))

    self._is_connected = False
    self._connected_subgraphs = []

    # If the given custom getter is a dictionary with a per-variable custom
    # getter, wrap it into a single custom getter.
    if isinstance(custom_getter, collections.Mapping):
      self._custom_getter = util.custom_getter_router(
          custom_getter_map=custom_getter,
          name_fn=lambda name: name[len(self.scope_name) + 1:])
    else:
      if not (custom_getter is None or callable(custom_getter)):
        raise TypeError("Given custom_getter is not callable.")
      self._custom_getter = custom_getter
    self._custom_getter = _maybe_wrap_custom_getter(
        _variable_tracking_custom_getter, self._custom_getter)

    self._template = tf.make_template(name,
                                      self._build_wrapper,
                                      create_scope_now_=True,
                                      custom_getter_=self._custom_getter)

    self._original_name = name
    self._unique_name = self._template.variable_scope.name.split("/")[-1]

    # Update __call__ and the object docstrings to enable better introspection.
    self.__doc__ = self._build.__doc__
    self.__call__.__func__.__doc__ = self._build.__doc__

    # Keep track of which graph this module has been connected to. Sonnet
    # modules cannot be connected to multiple graphs, as transparent variable
    # sharing is impossible in that case.
    self._graph = None

    # Container for all variables created in this module and its sub-modules.
    self._all_variables = set([])

  def _build_wrapper(self, *args, **kwargs):
    """Function which will be wrapped in a Template to do variable sharing.

    Passes through all arguments to the _build method, and returns the
    corresponding outputs, plus the name_scope generated by this call of the
    template.

    Args:
      *args: args list for self._build
      **kwargs: kwargs dict for self._build

    Returns:
      A tuple containing (output from _build, scope_name).
    """
    output = self._build(*args, **kwargs)
    # Make a dummy subscope to check the name scope we are in. We could read
    # the name scope from one of the outputs produced, except that the outputs
    # could have been produced from a subscope instantiated by the build
    # function, for example if inner modules are present. Calling name_scope
    # here and creating a new subscope guarantees we get the right answer.
    # Because we don't create an ops inside this dummy scope, no extra memory
    # will be consumed.
    with tf.name_scope("dummy") as scope_name:
      this_scope_name = scope_name[:-len("/dummy/")]
    return output, this_scope_name

  def _check_init_called(self):
    """Checks that the base class's __init__ method has been called.

    Raises:
      NotInitializedError: `AbstractModule.__init__` has not been called.
    """
    try:
      self._template
    except AttributeError:
      raise NotInitializedError("You may have forgotten to call super at the "
                                "start of %s.__init__."
                                % self.__class__.__name__)

  def _set_module_info(self):
    """Creates a `ModuleInfo` and adds it to the graph collections."""
    self._module_info = base_info.ModuleInfo(
        module_name=self.module_name,
        scope_name=self.scope_name,
        class_name="{}.{}".format(
            self.__class__.__module__, self.__class__.__name__),
        connected_subgraphs=self._connected_subgraphs)
    self._graph.add_to_collection(base_info.SONNET_COLLECTION_NAME,
                                  self._module_info)

  def _check_same_graph(self):
    """Checks that the module is not being connect to multiple Graphs.

    An instance of a Sonnet module 'owns' the variables it contains, and permits
    seamless variable sharing. As such, connecting a single module instance to
    multiple Graphs is not possible - this function will raise an error should
    that occur.

    Raises:
      DifferentGraphError: if the module is connected to a different Graph than
        it was previously used in.
    """
    current_graph = tf.get_default_graph()
    if self._graph is None:
      self._graph = current_graph
      self._set_module_info()
    elif self._graph != current_graph:
      raise DifferentGraphError("Cannot connect module to multiple Graphs.")

  @abc.abstractmethod
  def _build(self, *args, **kwargs):
    """Add elements to the Graph, computing output Tensors from input Tensors.

    Subclasses must implement this method, which will be wrapped in a Template.

    Args:
      *args: Input Tensors.
      **kwargs: Additional Python flags controlling connection.

    Returns:
      output Tensor(s).
    """

  @contextlib.contextmanager
  def _capture_variables(self):
    """Adds variables used by this module to self._all_variables.

    Upon entering this context manager the module adds itself onto the top
    of the module call stack. Any variables created with `tf.get_variable()`
    inside `_build()` or `_enter_variable_scope()` while this module is on top
    of the call stack will be added to `self._all_variables`.

    Before exiting the context the module removes itself from the top of the
    call stack, and adds all of the variables in `self._all_variables` to its
    parent module (the new top) of the call stack.

    Yields:
      Nothing, the yield just transfers focus back to the inner context.
    """
    module_stack = _MODULE_STACKS.setdefault(self._graph, [])
    module_stack.append(self)
    try:
      # In eager mode, the template store keeps references to created variables
      # such that they survive even if there are no references to them in
      # Python code. Variables added to an eager template store are also added
      # to TensorFlow global collections (unlike regular variables created in
      # eager mode).
      # Ideally move re-entering store into TF's tpl.variable_scope.
      if tf.executing_eagerly():
        with self._template._template_store.as_default():  # pylint:disable=protected-access
          yield
      else:
        yield
    finally:
      # Remove `self` from `module_stack`, this happens as part of cleanup
      # even if an error is raised.
      module_stack.pop()

    if module_stack:
      # Peek into the stack to add created variables to the parent
      parent_module = module_stack[-1]
      parent_module._all_variables.update(self._all_variables)  # pylint: disable=protected-access

  def _add_connected_subgraph(self, call_method, outputs, subgraph_name_scope,
                              *inputs_args, **inputs_kwargs):
    """Adds a newly connected subgraph.

    Args:
      call_method: the function used to connect this Sonnet module to the graph.
      outputs: `call_method` outputs.
      subgraph_name_scope: name scope of the newly connected subgraph.
      *inputs_args: `self._build` inputs `*args`.
      **inputs_kwargs: `self._build` inputs `*kwargs`.
    """
    build_inputs = inspect.getcallargs(call_method,
                                       *inputs_args, **inputs_kwargs)

    # "self" should normally be in `build_inputs` but some people are decorating
    # their `_build` function with `memoize`, in which case the function
    # signature doesn't contain `self` anymore.

    if "self" in build_inputs:
      del build_inputs["self"]

    connected_subgraph = base_info.ConnectedSubGraph(
        module=self, name_scope=subgraph_name_scope,
        inputs=build_inputs,
        outputs=outputs)
    self._connected_subgraphs.append(connected_subgraph)

    connection_observer_stack = _CONNECTION_OBSERVER_STACKS.setdefault(
        self._graph, [])
    for observer in connection_observer_stack:
      observer(connected_subgraph)

  def __call__(self, *args, **kwargs):
    """Operator overload for calling.

    This is the entry point when users connect a Module into the Graph. The
    underlying _build method will have been wrapped in a Template by the
    constructor, and we call this template with the provided inputs here.

    Args:
      *args: Arguments for underlying _build method.
      **kwargs: Keyword arguments for underlying _build method.

    Returns:
      The result of the underlying _build method.
    """
    self._check_init_called()
    self._check_same_graph()
    with self._capture_variables():
      outputs, subgraph_name_scope = self._template(*args, **kwargs)
    self._is_connected = True
    if not tf.executing_eagerly():
      # In eager mode the module is called a lot more frequently than in graph
      # mode (for each training step) and so we don't keep track of connected
      # subgraphs (since there will be orders of magnitude more of them).
      self._add_connected_subgraph(self._build, outputs, subgraph_name_scope,
                                   *args, **kwargs)
    return outputs

  @property
  def name_scopes(self):
    """Returns a tuple of all name_scopes generated by this module."""
    if tf.executing_eagerly():
      raise NotSupportedError(
          "The name_scopes property is not supported in eager mode.")
    return tuple(subgraph.name_scope for subgraph in self._connected_subgraphs)

  @property
  def variable_scope(self):
    """Returns the variable_scope declared by the module.

    It is valid for library users to access the internal templated
    variable_scope, but only makes sense to do so after connection. Therefore we
    raise an error here if the variable_scope is requested before connection.

    The only case where it does make sense to access the variable_scope before
    connection is to get the post-uniquification name, which we support using
    the separate .scope_name property.

    Returns:
      variable_scope: `tf.VariableScope` instance of the internal `tf.Template`.

    Raises:
      NotConnectedError: If the module is not connected to the Graph.
    """
    self._ensure_is_connected()
    return self._template.variable_scope

  @property
  def scope_name(self):
    """Returns the full name of the Module's variable scope."""
    return self._template.variable_scope.name

  @property
  def module_name(self):
    """Returns the name of the Module."""
    return self._unique_name

  @property
  def is_connected(self):
    """Returns true iff the Module been connected to the Graph at least once."""
    return self._is_connected

  @property
  def graph(self):
    """Returns the Graph instance which the module is connected to, or None."""
    return self._graph

  @property
  def connected_subgraphs(self):
    """Returns the subgraphs created by this module so far."""
    if tf.executing_eagerly():
      raise NotSupportedError(
          "Connected sub-graphs are not tracked in eager mode.")
    return tuple(self._connected_subgraphs)

  @property
  def last_connected_subgraph(self):
    """Returns the last subgraph created by this module.

    Returns:
      The last connected subgraph.

    Raises:
      NotConnectedError: If the module is not connected to the Graph.
    """
    if tf.executing_eagerly():
      raise NotSupportedError(
          "Connected sub-graphs are not tracked in eager mode.")
    self._ensure_is_connected()
    return self._connected_subgraphs[-1]

  @classmethod
  def get_possible_initializer_keys(cls):
    """Returns the keys the dictionary of variable initializers may contain.

    This provides the user with a way of knowing the initializer keys that are
    available without having to instantiate a sonnet module. Subclasses may
    override this class method if they need additional arguments to determine
    what initializer keys may be provided.

    Returns:
      Set with strings corresponding to the strings that may be passed to the
          constructor.
    """
    return getattr(cls, "POSSIBLE_INITIALIZER_KEYS", set())

  def _ensure_is_connected(self):
    """Raise an Error if the module has not been connected yet.

    Until the module is connected into the Graph, any variables created do
    not exist yet and cannot be created in advance due to not knowing the size
    of the input Tensor(s). This assertion ensures that any variables contained
    in this module must now exist.

    Raises:
      NotConnectedError: If the module is not connected to the Graph.
    """
    if not self.is_connected:
      raise NotConnectedError(
          "Variables in {} not instantiated yet, __call__ the module "
          "first.".format(self.scope_name))

  # pylint: disable=g-doc-return-or-yield
  @contextlib.contextmanager
  def _enter_variable_scope(self, reuse=None):
    """Returns a contextlib.contextmanager to enter the internal variable scope.

    This is useful for situations where submodules must be declared in the
    constructor, or somewhere else that is not called under the `_build` method.
    If such a case arises, calling `with self._enter_variable_scope():` will
    cause the variables in the submodule to be correctly scoped.

    An example justification for this is to allow the `Transposable` interface
    to be implemented - you might want to construct all the submodules at
    construction time so that you can call `.transpose()` and connect the
    result of that before connecting the non-transposed module.

    ```python
    class SomeModule(snt.AbstractModule):
      def __init__(self, name="some_module"):
        super(SomeModule, self).__init__(name=name)
        with self._enter_variable_scope():
          # We need to construct this submodule before we get to the _build
          # method, for some reason.
          self._sub_mod = snt.SomeSubmodule(name="some_submodule")

      def _build(self, input):
        # Connect to the already constructed submodule.
        return self._sub_mod(input)
    ```

    If you omit this then the submodule and parent module will appear to
    be "side by side" rather than nested when viewed in the Graph viewer, and
    functions such as `snt.get_variables_in_module()` or the `get_variables()`
    method will not know about variables defined in the submodule.

    Args:
      reuse: Boolean passed to `tf.variable_scope`.

    Yields:
      The variable_scope inside the template.
    """
    self._check_init_called()
    self._check_same_graph()
    with self._capture_variables():
      with tf.variable_scope(self._template.variable_scope, reuse=reuse) as vs:
        yield vs
  # pylint: enable=g-doc-return-or-yield

  def get_variables(self, collection=tf.GraphKeys.TRAINABLE_VARIABLES):
    """Returns tuple of `tf.Variable`s declared inside this module.

    Note that this operates by searching this module's variable scope,
    and so does not know about any modules that were constructed elsewhere but
    used inside this module.

    This method explicitly re-enters the Graph which this module has been
    connected to.

    Args:
      collection: Collection to restrict query to. By default this is
        `tf.Graphkeys.TRAINABLE_VARIABLES`, which doesn't include non-trainable
        variables such as moving averages.

    Returns:
      A tuple of `tf.Variable` objects.

    Raises:
      NotConnectedError: If the module is not connected to the Graph.
    """
    self._ensure_is_connected()
    # Explicitly re-enter Graph, in case the module is being queried with a
    # different default Graph from the one it was connected to. If this was not
    # here then querying the variables from a different graph scope would
    # produce an empty tuple.
    with self._graph.as_default():
      return util.get_variables_in_scope(
          self.variable_scope, collection=collection)

  def get_all_variables(self, collection=tf.GraphKeys.TRAINABLE_VARIABLES):
    """Returns all `tf.Variable`s used when the module is connected.

    See the documentation for `AbstractModule._capture_variables()` for more
    information.

    Args:
      collection: Collection to restrict query to. By default this is
        `tf.Graphkeys.TRAINABLE_VARIABLES`, which doesn't include non-trainable
        variables such as moving averages.

    Returns:
      A sorted (by variable name) tuple of `tf.Variable` objects.

    Raises:
      NotConnectedError: If the module is not connected to the Graph.
    """
    self._ensure_is_connected()
    collection_variables = set(tf.get_collection(collection))
    # Return variables in self._all_variables that are in `collection`
    return tuple(
        sorted(
            self._all_variables & collection_variables, key=lambda v: v.name))

  def __getstate__(self):
    raise NotSupportedError(
        "Sonnet AbstractModule instances cannot be serialized. You should "
        "instead serialize all necessary configuration which will allow "
        "modules to be rebuilt.")


@six.add_metaclass(abc.ABCMeta)
class Transposable(object):
  """Transposable module interface.

    The Transposable interface requires that transposable modules implement
    a method called `transpose`, returning a module that is the transposed
    version of the one the method is called on.
    Calling the method twice should return a module with the same specifications
    as the original module.

    When implementing a transposable module, special care is required to make
    sure that parameters needed to instantiate the module are provided as
    functions whose invocation is deferred to graph construction time.

    For example, in Linear we might want to call:

    ```python
    linear = snt.Linear(name="linear", output_size=output_size)
    linear_transpose = linear.transpose()
    ```

    where the output_size for linear_transpose is not known yet, as linear is
    not yet connected to the graph: output_size is passed to linear_transpose's
    constructor as a lambda returning linear.input_size. The lambda will return
    the correct value once linear is given an input.
    Notice that linear_transpose's output_size value does not need to be defined
    until the module is connected to the graph.
  """

  @abc.abstractmethod
  def transpose(self, name=None, **kwargs):
    """Builds and returns transposed version of module.

    Args:
      name: Name of the transposed module.
      **kwargs: Additional Python flags controlling transposition.

    Returns:
      Transposed version of the module.
    """

  @abc.abstractmethod
  def input_shape(self):
    """Returns shape of input `Tensor` passed at last call to `build`."""


class Module(AbstractModule):
  """Module wrapping a function provided by the user."""

  def __init__(self, build, custom_getter=None, name=None):
    """Constructs a module with a given build function.

    The Module class can be used to wrap a function assembling a network into a
    module.

    For example, the following code implements a simple one-hidden-layer MLP
    model by defining a function called make_model and using a Module instance
    to wrap it.

    ```python
    def make_model(inputs):
      lin1 = snt.Linear(name="lin1", output_size=10)(inputs)
      relu1 = tf.nn.relu(lin1, name="relu1")
      lin2 = snt.Linear(name="lin2", output_size=20)(relu1)
      return lin2

    model = snt.Module(name='simple_mlp', build=make_model)
    outputs = model(inputs)
    ```

    The `partial` package from `functools` can be used to bake configuration
    parameters into the function at construction time, as shown in the following
    example.

    ```python
    from functools import partial

    def make_model(inputs, output_sizes):
      lin1 = snt.Linear(name="lin1", output_size=output_sizes[0])(inputs)
      relu1 = tf.nn.relu(lin1, name="relu1")
      lin2 = snt.Linear(name="lin2", output_size=output_sizes[1])(relu1)
      return lin2

    model = snt.Module(name='simple_mlp',
                       build=partial(make_model, output_size=[10, 20])
    outputs = model(inputs)
    ```

    Args:
      build: Callable to be invoked when connecting the module to the graph.
          The `build` function is invoked when the module is called, and its
          role is to specify how to add elements to the Graph, and how to
          compute output Tensors from input Tensors.
          The `build` function signature can include the following parameters:
            *args - Input Tensors.
            **kwargs - Additional Python parameters controlling connection.
      custom_getter: Callable or dictionary of callables to use as
          custom getters inside the module. If a dictionary, the keys
          correspond to regexes to match variable names. See the
          `tf.get_variable` documentation for information about the
          custom_getter API.
      name: Module name. If set to `None` (the default), the name will be set to
          that of the `build` callable converted to `snake_case`. If `build` has
          no name, the name will be 'module'.

    Raises:
      TypeError: If build is not callable.
      TypeError: If a given `custom_getter` is not callable.
    """
    if not callable(build):
      raise TypeError("Input 'build' must be callable.")
    if name is None:
      name = util.name_for_callable(build)
    super(Module, self).__init__(custom_getter=custom_getter, name=name)
    self._build_function = build

  def _build(self, *args, **kwargs):
    """Forwards call to the passed-in build function."""
    return self._build_function(*args, **kwargs)
