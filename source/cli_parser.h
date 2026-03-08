#ifndef CLI_PARSER_H
#define CLI_PARSER_H

#include <string>

enum class RunMode { UNKNOWN, SIM_RECORD_FIELD, SIM_RECORD_RESPONSE, SIM_VIZ, VIZ_RECORD };

struct CliArgs {
    RunMode mode = RunMode::UNKNOWN;
    std::string experiment_name = "";
    std::string config_path = "";
    std::string playback_file = "";
    int playback_delay = 33;
};

CliArgs parse_cli_args(int argc, char* argv[]);

#endif // CLI_PARSER_H
