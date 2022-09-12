/*******************************************************************************
* Copyright 2020 Intel Corporation.
*
* This software and the related documents are Intel copyrighted  materials,  and
* your use of  them is  governed by the  express license  under which  they were
* provided to you (License).  Unless the License provides otherwise, you may not
* use, modify, copy, publish, distribute,  disclose or transmit this software or
* the related documents without Intel's prior written permission.
*
* This software and the related documents  are provided as  is,  with no express
* or implied  warranties,  other  than those  that are  expressly stated  in the
* License.
*******************************************************************************/

#ifndef __SYCL_DEVICE_COMMON_FOR_EXAMPLES_HPP__
#define __SYCL_DEVICE_COMMON_FOR_EXAMPLES_HPP__

#include <list>
#include <map>

#include <CL/sycl.hpp>

//
// custom helpers for referring to different sycl devices
//
enum my_sycl_device_types {
    host_device,
    cpu_device,
    gpu_device
};

std::map<my_sycl_device_types, std::string>
sycl_device_names = { {host_device, "Host"},
                      {cpu_device, "CPU"},
                      {gpu_device, "GPU"}  };

//
// users can set environment flags like SYCL_DEVICES_{all,host,cpu,gpu} to specify
// which devices to run on
//
void set_list_of_devices(std::list<my_sycl_device_types> & list_of_devices) {
#if defined(SYCL_DEVICES_all)
    list_of_devices.push_back(host_device);
    list_of_devices.push_back(cpu_device);
    list_of_devices.push_back(gpu_device);
#else
#if defined(SYCL_DEVICES_host) 
    list_of_devices.push_back(host_device);
#endif
#if defined(SYCL_DEVICES_cpu)
    list_of_devices.push_back(cpu_device);
#endif
#if defined(SYCL_DEVICES_gpu)
    list_of_devices.push_back(gpu_device);
#endif
#endif
}

//
// sets device using *_selector() functionality and returns whether it was successful at finding it
//
void get_sycl_device(sycl::device &my_dev, bool & my_dev_is_found, my_sycl_device_types & desired_sycl_device) {

    my_dev_is_found = true;

    try {
        switch (desired_sycl_device) {
            case host_device:
                my_dev = sycl::device(sycl::host_selector());
                break;
            case cpu_device:
                my_dev = sycl::device(sycl::cpu_selector());
                break;
            case gpu_device:
                my_dev = sycl::device(sycl::gpu_selector());
                break;
        }
    } catch (...) {
        my_dev_is_found = false;
    }
}

bool isDoubleSupported(sycl::device my_dev) {
    return my_dev.get_info<sycl::info::device::double_fp_config>().size() != 0;
}

#endif // __SYCL_DEVICE_COMMON_FOR_EXAMPLES_HPP__
