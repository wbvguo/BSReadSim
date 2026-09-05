#ifndef OPTION_H
#define OPTION_H

// Command-line options, "--name value" (booleans are true or false). The Python package
// passes every option (the names follow its settings); run on its own, htsim needs only
// --reference and uses the defaults of the parameter structs for the rest. parse_options
// also checks how the options combine, for `command` (0 for a simulation run).

#include "struct.h"

// parse argv (after the program or the subcommand name) into the parameters
void parse_options(int argc, char *argv[], const char *command,
                   expt_param *expt_set, mut_param *mut_set, meth_param *meth_set);

#endif
