#!/bin/bash
# Script for to get machine information before doing the experiment

set +e # Don't fail fast since some information is maybe not available

title="Experiment results"
inputfile=""
host="$(hostname | sed 's/[0-9]*//g' | cut -d'.' -f1)"
help_script()
{
      cat << EOF
Usage: $0 [options] outputfile.org

Script for to get machine information before doing the experiment

OPTIONS:
   -h      Show this message
   -t      Title of the output file
EOF
}
# Parsing options
while getopts "t:s:i:h" opt; do
      case $opt in
          t)
                  title="$OPTARG"
                        ;;
                          h)
                                  help_script
                                        exit 4
                                              ;;
                                                \?)
                                                        echo "Invalid option: -$OPTARG"
                                                              help_script
                                                                    exit 3
                                                                          ;;
                                                                              esac
                                                                            done

                                                                            shift $((OPTIND - 1))
                                                                            filedat=$1
                                                                            if [[ $# != 1 ]]; then
                                                                                  echo 'ERROR!'
                                                                                      help_script
                                                                                          exit 2
                                                                            fi

                                                                            ##################################################
                                                                            # Preambule of the output file
                                                                            echo "#+TITLE: $title" >> $filedat
                                                                            echo "#+DATE: $(eval date)" >> $filedat
                                                                            echo "#+AUTHOR: $(eval whoami)" >> $filedat
                                                                            echo "#+MACHINE: $(eval hostname)" >> $filedat
                                                                            echo "#+FILE: $(eval basename $filedat)" >> $filedat
                                                                            if [[ -n "$inputfile" ]]; 
                                                                            then
                                                                                  echo "#+INPUTFILE: $inputfile" >> $filedat
                                                                            fi
                                                                            echo " " >> $filedat 

                                                                            ##################################################
                                                                            # Collecting metadata
                                                                            echo "* MACHINE INFO:" >> $filedat

                                                                            echo "** PEOPLE LOGGED WHEN EXPERIMENT STARTED:" >> $filedat
                                                                            who >> $filedat
                                                                            echo "############################################" >> $filedat

                                                                            echo "** ENVIRONMENT VARIABLES:" >> $filedat
                                                                            env >> $filedat
                                                                            echo "############################################" >> $filedat

                                                                            echo "** HOSTNAME:" >> $filedat
                                                                            hostname >> $filedat
                                                                            echo "############################################" >> $filedat

                                                                            if [[ -n $(command -v lstopo) ]];
                                                                            then
                                                                                  echo "** MEMORY HIERARCHY:" >> $filedat
                                                                                      lstopo --of console >> $filedat
                                                                                          echo "############################################" >> $filedat
                                                                            fi

                                                                            if [ -f /proc/cpuinfo ];
                                                                            then
                                                                                  echo "** CPU INFO:" >> $filedat
                                                                                      cat /proc/cpuinfo >> $filedat
                                                                                          echo "############################################" >> $filedat
                                                                            fi

                                                                            if [ -f /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor ];
                                                                            then
                                                                                  echo "** CPU GOVERNOR:" >> $filedat
                                                                                      ONLINECPUS=$(for CPU in $(find /sys/devices/system/cpu/ | grep cpu[0-9]*$); do [[ $(cat $CPU/online) -eq 1 ]] && echo $CPU; done | grep cpu[0-9]*$ | sed 's/.*cpu//')
                                                                                          for PU in ${ONLINECPUS}; do
                                                                                                   echo -n "CPU frequency for cpu${PU}: " >> $filedat
                                                                                                          cat /sys/devices/system/cpu/cpu${PU}/cpufreq/scaling_governor >> $filedat
                                                                                                              done
                                                                                                                  echo "############################################" >> $filedat
                                                                            fi

                                                                            if [ -f /sys/devices/system/cpu/cpu0/cpufreq/scaling_cur_freq ];
                                                                            then
                                                                                  echo "** CPU FREQUENCY:" >> $filedat
                                                                                      ONLINECPUS=$(for CPU in $(find /sys/devices/system/cpu/ | grep cpu[0-9]*$); do [[ $(cat $CPU/online) -eq 1 ]] && echo $CPU; done | grep cpu[0-9]*$ | sed 's/.*cpu//')
                                                                                          for PU in ${ONLINECPUS}; do
                                                                                                   echo -n "CPU frequency for cpu${PU}: " >> $filedat
                                                                                                         cat /sys/devices/system/cpu/cpu${PU}/cpufreq/scaling_cur_freq >> $filedat
                                                                                                             done
                                                                                                                 echo "############################################" >> $filedat
                                                                            fi

                                                                            if [ -f /usr/bin/cpufreq-info ];
                                                                            then
                                                                                  echo "** CPUFREQ_INFO" >> $filedat
                                                                                      cpufreq-info >> $filedat
                                                                                          echo "############################################" >> $filedat
                                                                            fi

                                                                            if [ -f /usr/bin/lspci ];
                                                                            then
                                                                                  echo "** LSPCI" >> $filedat
                                                                                      lspci >> $filedat
                                                                                          echo "############################################" >> $filedat
                                                                            fi

                                                                            if [ -f /usr/bin/ompi_info ];
                                                                            then
                                                                                  echo "** OMPI_INFO" >> $filedat
                                                                                      ompi_info --all >> $filedat
                                                                                          echo "############################################" >> $filedat
                                                                            fi

                                                                            if [ -f /sbin/ifconfig ];
                                                                            then
                                                                                  echo "** IFCONFIG" >> $filedat
                                                                                      /sbin/ifconfig >> $filedat
                                                                                          echo "############################################" >> $filedat
                                                                            fi

                                                                            if [[ -n $(command -v nvidia-smi) ]];
                                                                            then
                                                                                  echo "** GPU INFO FROM NVIDIA-SMI:" >> $filedat
                                                                                      nvidia-smi -q >> $filedat
                                                                                          echo "############################################" >> $filedat
                                                                            fi 

                                                                            if [ -f /proc/version ];
                                                                            then
                                                                                  echo "** LINUX AND GCC VERSIONS:" >> $filedat
                                                                                      cat /proc/version >> $filedat
                                                                                          echo "############################################" >> $filedat
                                                                            fi

                                                                            if [[ -n $(command -v module) ]];
                                                                            then
                                                                                  echo "** MODULES:" >> $filedat
                                                                                      module list 2>> $filedat
                                                                                          echo "############################################" >> $filedat
                                                                            fi

                                                                            echo "** TCP PARAMETERS" >> $filedat
                                                                            FILES="/proc/sys/net/core/rmem_max \
                                                                              /proc/sys/net/core/wmem_max \
                                                                              /proc/sys/net/core/rmem_default \
                                                                              /proc/sys/net/core/wmem_default \
                                                                              /proc/sys/net/core/netdev_max_backlog \
                                                                              /proc/sys/net/ipv4/tcp_rmem \
                                                                              /proc/sys/net/ipv4/tcp_wmem \
                                                                              /proc/sys/net/ipv4/tcp_mem"

                                                                            for FILE in $FILES; do
                                                                                  echo "cat $FILE"
                                                                                      cat $FILE
                                                                                    done >> $filedat
