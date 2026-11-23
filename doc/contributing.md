title: Contributing to AGILE and its documentation

# Introduction {: #contrib-intro }
This page describes how you can contribute code and documentation to AGILE.

# Testing your changes {: #contrib-testing }
Having made some changes, the first thing you can check is whether the existing 
test suite still works.   This is done as
follows:

```text
cd tests
make all
```

This will run the prepared test cases and compare their output to previously
stored results, and report errors when differences are larger than some
threshold. Alternatively, you can run tests for individual physics modules:

```text
cd tests
make hd
make mhd
...
```

## Adding tests {: #contrib-adding-tests }
To add new tests, you need:

* A `mod_usr.fpp`
* A suitable `.par` file so that your simulation runs quickly (in say 1 to 10 seconds)
* A file `test.make` following the template below

These files should be placed in the corresponding physics folder, e.g. `tests/hd`
for hydro problems. Assuming your `.par` file is called `test.par`, it should contain the following options in its filelist (see [Setting parameters](par.html)):

```fortran
    base_filename='my_test'
    typefilelog='regression_test'
```

You can then use the following template for the `test.make` file:

```make
    # Template for test.make files

    TESTS := test.log
    
    include ../../test_rules.make
    
    # Generate dependency rules for the tests
    test.log: test.par
```

Of course, you can change `test.log` in the above files to something more
meaningful, e.g. `hd_dust`. Now you can compile and run your tests with:

```text
make -f test.make
```

They should all fail the first time, but there should now be log files of the
"regression_test" type, which contain the volume averages of the variables and
the variables squared over time. If you are happy with the results, store them
in a folder called `correct_output`:

```text
mkdir correct_output
cp *.log correct_output/
```

Now `make -f test.make` should pass. You can then add your test folder to the
corresponding directory in `tests/Makefile` file, so that the test is
automatically performed. Finally, commit and push your changes.

# Style guide {: #contrib-style }
To steadily improve AGILE, it would be good if new contributions take into
account our [style guide](code_style_guide.html). Most importantly:

* Choose meaningful names for variables, functions etc.
* For non-trivial blocks of code or routines: describe what they are supposed to do

# Working with git {: #contrib-git }
Git is a version control system that you can use for code and other 'line-based'
documents, such as LaTeX files. There are many tutorials on git, so it is
probably best to look at a couple of them and pick one you like, for example:

* https://help.github.com/articles/good-resources-for-learning-git-and-github/
* https://www.atlassian.com/git/tutorials/
* https://git-scm.com/docs/gittutorial
* http://rogerdudler.github.io/git-guide/

You can ask for access to the AGILE Github repository so that you can
`push` your changes to it. You can also clone the repository and/or file a merge
request.

# Contributing documentation {: #contrib-doc }
The [documentation](documentation.html) page explains how to write documentation
for AGILE.
