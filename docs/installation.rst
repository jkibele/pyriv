.. highlight:: shell

============
Installation
============
** Note: Python 2.6 is not supported due to a bug in 2.6 that prevents the use of partials. In the future, the reliance on partials for parallelization should be replaced with a wrapper function or something of the sort.**

Stable release
--------------

To install pyriv, run this command in your terminal:

.. code-block:: console

    $ pip install pyriv

This is the preferred method to install pyriv, as it will always install the most recent stable release. 

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/


From sources
------------

The sources for pyriv can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/jkibele/pyriv

Or download the `tarball`_:

.. code-block:: console

    $ curl  -OL https://github.com/jkibele/pyriv/tarball/master

Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ python setup.py install


.. _Github repo: https://github.com/jkibele/pyriv
.. _tarball: https://github.com/jkibele/pyriv/tarball/master
