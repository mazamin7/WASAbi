#include "cli_parser.h"

using namespace std;

CliArgs parse_cli_args(int argc, char* argv[]) {
    CliArgs args;
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--mode" && i + 1 < argc) {
            string m = argv[++i];
            if (m == "sim-record-field") args.mode = RunMode::SIM_RECORD_FIELD;
            else if (m == "sim-record-response") args.mode = RunMode::SIM_RECORD_RESPONSE;
            else if (m == "sim-viz") args.mode = RunMode::SIM_VIZ;
            else if (m == "viz-record") args.mode = RunMode::VIZ_RECORD;
        }
        else if (arg == "--config" && i + 1 < argc) args.config_path = argv[++i];
        else if (arg == "--playback-file" && i + 1 < argc) args.playback_file = argv[++i];
        else if (arg == "--playback-delay" && i + 1 < argc) args.playback_delay = atoi(argv[++i]);
    }
    return args;
}
