#!/usr/bin/python3.5
"""
Checks whether another module, class or funcion
is documented (has a docstring).
"""
import os
import sys
import inspect

from types import ModuleType, MethodType, FunctionType
from operator import itemgetter

import simos


class Check:
  """
  Check whether `obj` and all its sub-objects are documented.
  Only documentable objects are checked.

  `__new__` returns a list of tuples (obj, has_doc).
  """

  def __new__(cls, obj):
    """
    Check whether `obj` and all its sub-objects are documented.
    Only documentable objects are checked.

    :param obj: object to be checked
    :return: a list of tuples (obj, has_doc)
    """
    self = object.__new__(cls)
    self.objs = []
    self.values = []
    self.check(obj)
    return list(zip(self.objs, self.values))


  def check(self, obj):
    """
    Check whether `obj` and all its sub-objects are documented.
    Only documentable objects are checked.
    ! Only updates `self.objs` and `self.values`;
    ! doesn't return anything.

    :param obj: object to be checked
    """

    obj = self.documentable(obj)

    if obj is None or obj in self.objs:
      return

    self.objs.append(obj)
    self.values.append(self.has_doc(obj))

    self.check_children(obj)


  @staticmethod
  def documentable(obj):
    """
    Return a documentable object from `obj`
    based on the type of `obj`. If type is:
     - module, type, or function => `obj` itself
     - classmethod or staticmethod => `obj.__func__`
     - else `None` (meaning the object is not documetable)

    :param obj: object to be converted toa  documentable object
    :return: documentable object or `None`
    """
    if isinstance(obj, (classmethod, staticmethod)):
      return obj.__func__

    if isinstance(obj, ModuleType) or \
       isinstance(obj, type) or \
       isinstance(obj, MethodType) or \
       (isinstance(obj, FunctionType)
        and obj.__name__ !=  '<lambda>'):
          return obj


  @staticmethod
  def has_doc(obj):
    """
    Check whether a documentable object is documented.

    :param obj: object to be checked
    :return: whether `obj` is documented (bool)
    """
    return obj.__doc__ is not None


  def check_children(self, obj):
    """
    If `obj` is a module, type, or function, this function checks
    whether all documentable sub-objects in `obj` are documented.

    :param obj: object, sub-objects of which are to be checked
    """
    if isinstance(obj, ModuleType) or\
       isinstance(obj, type):
         self.check_by_dict(obj)

  def check_by_dict(self, obj):
    """
    Check whether all documentable objects in `obj.__dict__`
    are documented.

    :param obj: object, sub-objects of which are to be checked
    """
    for k,v in obj.__dict__.items():
      self.check(v)


class anykey_defaultdict:
  """
  `defaultdict` supporting any values as keys.
  """

  def __init__(self, factory):
    """
    Init self.

    :param factory: factory for default values
    """
    self.factory = factory
    self.keys   = []
    self.values = []

  def __getitem__(self, key):
    """
    Return item for `key`. If not present, use `self.factory`
    to create the item, save the item in self and return it.
    """
    if key not in self.keys:
      self.keys.append(key)
      value = self.factory()
      self.values.append(value)
    else:
      index = self.keys.index(key)
      value = self.values[index]
    return value

  def __setitem__(self, key, value):
    """
    Assign `value` to `key`.
    """
    if key not in self.keys:
      self.keys.append(key)
      self.values.append(value)
    else:
      index = self.keys.index(key)
      self.values[index] = value

  def __delitem__(self, key):
    """
    Remove item corresponding to `key` from self.

    The order of items is not preserved after deletion,
    so indices obtained before deleting an item
    might be incorrect after the deletion.
    """
    if key not in self.keys:
      raise KeyError(key)

    index = self.key.index(key)

    self.keys[index] = self.keys[-1]
    del self.keys[-1]

    self.values[index] = self.values[-1]
    del self.values[-1]


  def set(self, index, value):
    """
    Assign `value` to item at `index`.

    If an item was deleted between obtaining the index
    and calling this method, the index might be incorrect.

    :param index: index at which the value is to be assigned
    :param value: value to be assigned
    """
    self.values[index] = value

  def __iter__(self):
    """
    Return iterator of values.
    """
    return iter(self.values)


def group_by(objs, *groupers):
  """
  Group `objs` by `groupers`.
  """
  if not groupers:
    return objs

  grouper = groupers[0]
  grouped = anykey_defaultdict(list)
  for obj in objs:
    group = grouper(obj)
    grouped[group].append(obj)

  for i,group in enumerate(grouped):
    grouped.set(i, group_by(group, *groupers[1:]))

  return list(zip(grouped.keys, grouped.values))


def file_grouper(obj):
  """
  Group by `inspect.getfile(obj)`.
  """
  try:
    return inspect.getfile(obj)
  except TypeError:
    return ''


def check_groupby_file_type(obj):
  """
  Check `obj`, group by file, select only undocumented objects.

  :param obj: object to check
  :return: list of undocumented objects grouped by (filepath, type)
  """
  obj_dpath = os.path.dirname(inspect.getfile(obj))

  def doesnt_have_doc(obj_hd):
    """
    Check if the object does not have a docstring.

    :param obj_hd: tuple containing the object and its documentation status
    :return: True if the object does not have a docstring, False otherwise
    """
    return not obj_hd[1]
  objs_without_doc = map(itemgetter(0), filter(doesnt_have_doc, Check(obj)))

  objs_grouped_by_file_type = group_by(objs_without_doc, file_grouper, type)

  def is_local(file_value):
    """
    Check if the file path starts with the object's directory path.

    :param file_value: tuple containing the file path and the value
    :return: True if the file path is local, False otherwise
    """
    return file_value[0].startswith(obj_dpath)

  target_objs = filter(is_local, objs_grouped_by_file_type)

  return list(target_objs)


def namegetter(obj):
  """
  Get name (string representation) of `obj`, as "qualified" as possible:
  - `obj.__qualname__` if available,
  - else `obj.__name__` if available,
  - else `str(obj)`.

  :param obj: object, name of which to return
  :return: name of `obj`
  """
  return getattr(obj, '__qualname__', None) \
      or getattr(obj, '__name__', None) \
      or str(obj)


def check(obj):
  """
  !IMPURE! Check whether `obj` and sub-objects are documented.

  Print undocumented objects from `obj`'s directory grouped by filepath
  and type (function, class, module,...), sorted by name.
  Return `True` if all objects are documented, `False` otherwise.

  :param obj: object to check
  :return: whether the object and all its sub-objects are documented
  """
  without_doc = check_groupby_file_type(obj)

  for fpath, list_kind_objs in without_doc:
    print('\x1b[1;31m%s\x1b[m' % fpath)

    for kind, objs in list_kind_objs:
      print('\t\x1b[1m%s\x1b[m' % kind.__qualname__)

      for obj_qualname in sorted(map(namegetter, objs)):
        print('\t\t%s' % obj_qualname)

  return not without_doc


def main(*args):
  """
  !IMPURE! Check all modules given in `args`.

  Print names of undocumented modules, classes and functions
  grouped by filepath and type, sorted by name.

  :param args: paths (Python or OS) to modules to be checked
  :return: number of modules which lack documentation
  """
  import os
  import importlib
  fails = 0

  for mod_fpath in sys.argv[1:]:
    mod_dpath = os.path.dirname(mod_fpath)

    if mod_fpath.endswith('.py'):
      mod_fpath = mod_fpath[:-len('.py')]

    mod_fname = os.path.basename(mod_fpath)

    sys.path.insert(0, mod_dpath)
    #print('\x1b[1m%s\x1b[m' % mod_pypath)
    fails += not check(importlib.import_module(mod_fname))
    del sys.path[0]

  return fails


print(check(simos))