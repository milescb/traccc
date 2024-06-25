/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// io
#include "traccc/io/read_cells.hpp"
#include "traccc/io/read_digitization_config.hpp"
#include "traccc/io/read_geometry.hpp"
#include "traccc/io/utils.hpp"
#include "traccc/io/write.hpp"

// algorithms
#include "traccc/ambiguity_resolution/greedy_ambiguity_resolution_algorithm.hpp"
#include "traccc/clusterization/clusterization_algorithm.hpp"

// performance
#include "traccc/performance/timer.hpp"

// options
#include "traccc/options/clusterization.hpp"
#include "traccc/options/detector.hpp"
#include "traccc/options/input_data.hpp"
#include "traccc/options/output_data.hpp"
#include "traccc/options/performance.hpp"
#include "traccc/options/program_options.hpp"

// Detray include(s).
#include "detray/core/detector.hpp"
// #include "detray/detectors/bfield.hpp"
#include "detray/io/frontend/detector_reader.hpp"
#include "detray/navigation/navigator.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s).
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <map>
#include <memory>

int seq_run(const traccc::opts::input_data& input_opts,
            const traccc::opts::detector& detector_opts,
            const traccc::opts::clusterization&) {

    // Memory resource used by the application.
    vecmem::host_memory_resource host_mr;

    // Read in the geometry.
    auto [surface_transforms, barcode_map] = traccc::io::read_geometry(
        detector_opts.detector_file,
        (detector_opts.use_detray_detector ? traccc::data_format::json
                                           : traccc::data_format::csv));

    using detector_type = detray::detector<detray::default_metadata,
                                           detray::host_container_types>;
    detector_type detector{host_mr};
    if (detector_opts.use_detray_detector) {
        // Set up the detector reader configuration.
        detray::io::detector_reader_config cfg;
        cfg.add_file(traccc::io::data_directory() +
                     detector_opts.detector_file);
        if (detector_opts.material_file.empty() == false) {
            cfg.add_file(traccc::io::data_directory() +
                         detector_opts.material_file);
        }
        if (detector_opts.grid_file.empty() == false) {
            cfg.add_file(traccc::io::data_directory() +
                         detector_opts.grid_file);
        }

        // Read the detector.
        auto det = detray::io::read_detector<detector_type>(host_mr, cfg);
        detector = std::move(det.first);
    }

    // Read the digitization configuration file
    auto digi_cfg =
        traccc::io::read_digitization_config(detector_opts.digitization_file);

    // Output stats
    uint64_t n_cells = 0;
    uint64_t n_modules = 0;
    uint64_t n_measurements = 0;

    // Algorithms
    traccc::host::clusterization_algorithm ca(host_mr);

    // Timers
    traccc::performance::timing_info elapsedTimes;

    // Loop over events
    for (unsigned int event = input_opts.skip;
         event < input_opts.events + input_opts.skip; ++event) {

        traccc::host::clusterization_algorithm::output_type
            measurements_per_event{&host_mr};

        {  // Start measuring wall time.
            traccc::performance::timer timer_wall{"Wall time", elapsedTimes};

            traccc::io::cell_reader_output readOut(&host_mr);

            {
                traccc::performance::timer timer{"Read cells", elapsedTimes};
                // Read the cells from the relevant event file
                traccc::io::read_cells(readOut, event, input_opts.directory,
                                       input_opts.format, &surface_transforms,
                                       &digi_cfg, barcode_map.get());
            }
            traccc::cell_collection_types::host& cells_per_event =
                readOut.cells;
            traccc::cell_module_collection_types::host& modules_per_event =
                readOut.modules;

            /*-------------------
                Clusterization
              -------------------*/

            {
                traccc::performance::timer timer{"Clusterization",
                                                 elapsedTimes};
                measurements_per_event =
                    ca(vecmem::get_data(cells_per_event),
                       vecmem::get_data(modules_per_event));

                auto measurements_size = measurements_per_event.size();
                std::cout << "Number of measurements: " << measurements_size << std::endl;

                for (std::size_t i = 0; i < 10; ++i) {
                    auto measurement = measurements_per_event.at(i);
                    std::cout << "Measurement ID: " << measurement.measurement_id << std::endl;
                    std::cout << "Local coordinates: [" << measurement.local[0] << ", " << measurement.local[1] << "]" << std::endl; 
                }
            }

            /*----------------------------
              Statistics
              ----------------------------*/

            n_modules += modules_per_event.size();
            n_cells += cells_per_event.size();
            n_measurements += measurements_per_event.size();

        }  // Stop measuring Wall time.

    }

    std::cout << "==> Statistics ... " << std::endl;
    std::cout << "- read     " << n_cells << " cells from " << n_modules
              << " modules" << std::endl;
    std::cout << "- created  " << n_measurements << " measurements. "
              << std::endl;
    std::cout << "==> Elapsed times...\n" << elapsedTimes << std::endl;

    return EXIT_SUCCESS;
}

// The main routine
//
int main(int argc, char* argv[]) {

    // Program options.
    traccc::opts::detector detector_opts;
    traccc::opts::input_data input_opts;
    traccc::opts::clusterization clusterization_opts;
    traccc::opts::performance performance_opts;
    traccc::opts::program_options program_opts{
        "Clusterization on the Host",
        {detector_opts, input_opts, clusterization_opts, performance_opts
        },
        argc,
        argv};

    // Run the application.
    return seq_run(input_opts, detector_opts, 
            clusterization_opts);
}