#!/bin/bash
sh -c "source ~/.bashrc.custom"
hadoop datanode &
hadoop namenode &
hadoop jobtracker &
hadoop tasktracker &
