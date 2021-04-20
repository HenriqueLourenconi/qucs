--
-- README
--
-- Copyright (C) 2003, 2005, 2011 Stefan Jahn <stefan@lkcc.org>
--               2020 Felix Salfelder
--
-- This is free software; you can redistribute it and/or modify
-- it under the terms of the GNU General Public License as published by
-- the Free Software Foundation; either version 3, or (at your option)
-- any later version.
--
-- This software is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
-- GNU General Public License for more details.
--
-- You should have received a copy of the GNU General Public License
-- along with this package; see the file COPYING.  If not, write to
-- the Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
-- Boston, MA 02110-1301, USA.
--

## WIP WIP WIP

This is the Qucs refactoring repository which will bring Qt5 support among
other things. Help make it happen with your donation.
[![donate](https://liberapay.com/assets/widgets/donate.svg "donate through lp")](https://liberapay.com/Gnucap_and_Qucs/donate)


## Description

Qucs was a circuit simulator with a graphical user interface. "Qucs, the GUI",
this package, was the graphical user interface for circuit simulation. The
software aimed to support various kinds of circuit simulation, e.g. DC, AC,
S-parameter and harmonic balance analysis.

The current work aims at enhancing the extensibility and maintainability of
"Qucs the GUI" and its user contributions. It consists of a library
implementing circuit related data structures and and a GUI for interactive
schematic editing and simulation. User data, such as circuit components,
specific schematic/netlist languages, visualisation tools as well as
simulator drivers can be loaded at runtime.

## Requirements

Qucs needs Qt version 5 or later. For simplicity, the Qt configuration is
currently obtained through pkg-config, which you need to install in addition
to Qt.

On Debian, install

 pkg-config
 build-essential
 qtscript5-dev
 libqt5svg5-dev

For development purposes, autotools are required to bootstrap, see below.

## Installation

OPTION I: Unpack the distribution tarball.

    $ tar xvzf qucs-<version>.tar.gz               (using GNU tar)
    $ gzip -cd qucs-<version>.tar.gz | tar xvf -   (using another tar)
    $ cd qucs-<version>

OPTION II: Prepare source from git

    $ git clone $repo_url # github or so.
    $ cd qucs
    $ ./bootstrap
    $ [..]

Configure and build binaries for your system.

    $ mkdir build
    $ cd build
    $ ../configure
    $ make

Run tests (optional).

    $ make check

Install Qucs (optional).

    $ make install

You need write permissions to the prefix of the installation (default:
/usr/local). See INSTALL for further details.

## use CMake instead

$ ln -s cmake.stale cmake # for now

$ mkdir cmake/build
$ cd cmake/build
$ cmake ..
$ make install # this needs work

or

$ mkdir cb
$ cd cb
$ cmake ../cmake
$ make install # this needs work

## Getting the latest Git snapshot

You can get the latest Qucs version from our Git repository. Please use an
official release if you want to work with Qucs. The Git version might not
always compile or run smoothly. Currently, the main git repository is hosted
on github.

## Debugging

Configure a debug build as in

$ ../configure --enable-debug

It is not normally recommended to install during development. Running the
programs directly from the build directory is possible, with some caveats.
Debugging is easiest with libtool, but the environment variables that are not
handled by autotools need to be set first.

$ cd $build_directory
$ . main/qucs       # set environment variables for running from builddir
$ libtool --mode=execute gdb --args main/qucs.real -i examples/bandpass.sch
