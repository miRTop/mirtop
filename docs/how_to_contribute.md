# Examples of contributions

## How to add a new sub-command

**You need first to clone and install the tool in [develop mode](installation.html)**

Let's say that you want to add a new operation to `mirtop`, for instance, similar to the `stats` command to work with sGFF3 files. Assume a `test` function for this exmaple to just read the file and print `Hello GFF3.`

* Create the folder inside `mirtop/test`. The create to empty files named:
 * `test.py`
 * `__init__.py`

* Modify the `test.py` file with this content:

```
from mirtop.gff.body import read_gff_line

import mirtop.libs.logger as mylog
logger = mylog.getLogger(__name__)

def test(args):
	for fn in args.files:
		_test(fn)
		logger.info("Hello GFF3: %s" % fn)


def _test(fn):
	logger.debug("I am going to read this file: %s" % fn)
	for line in fn:
		read_gff_line(line)
```

* Choose a sub_command name, in this case: `test`.
* Add the arguments function at the end of this file: https://github.com/miRTop/mirtop/blob/dev/mirtop/libs/parse.py, using a naming following `add_subparser_test`. 

```
def add_subparser_test(subparsers):
    parser = subparsers.add_parser("test", help="test function")
    parser.add_argument("files", nargs="*", help="GFF/GTF files.")
    parser = _add_debug_option(parser)
    return parser
```

* Add the function name to `parse_cl` function, at the end of the `sub_cmds` array.

```
    sub_cmds = {"gff": add_subparser_gff,
                "stats": add_subparser_stats,
                "compare": add_subparser_compare,
                "target": add_subparser_target,
                "simulator": add_subparser_simulator,
                "counts": add_subparser_counts,
                "export": add_subparser_export,
                "test": add_subparser_test
                }
```

* To get the function re-directed from the command line when you use the sub_cmd name, add a line to the `command_line.py` file, adding another `else` statement:

```
	elif "test" in kwargs:
        logger.info("Run test.")
        test(kwargs["args"])
```

* The function you use to link to the operation added need to be imported at the beginning. Let's say that the `test` function is at `mirtop/test/test.py`:

```
from mirtop.test import test
```

Try the new operation:

```
mirtop test data/examples/correct_file.gff
```

## Add a unit test

## for the internal function


Add to the end of `test/test_functions.py`, but inside `class FunctionsTest(unittest.TestCase):` this code:


```
    @attr(fn_test=True)
    def test_function_test(self):
        from mirtop import test
        test._test("data/examples/gff/correct_file.gff")

```

## for the sub-command

Add to the end of `test/test_function.py`, but inside `class AutomatedAnalysisTest(unittest.TestCase):` this code:

```
    @attr(cmd_test=True)
    def test_srnaseq_annotation_bam(self):
        """Run test analysis
        """
        with make_workdir():
            clcode = ["mirtop",
                      "test",
                      "../../data/examples/gff/correct_file.gff"]
            print ""
            print " ".join(clcode)
            subprocess.check_call(clcode)

```

## test the unit

**nose is needed: pip install nose**

Run the function test from the top parent folder:

```
./run_test.sh fn_test
```

Run the command test from the top parent folder:

```
./run_test.sh cmd_test
```