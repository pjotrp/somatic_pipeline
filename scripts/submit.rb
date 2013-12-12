#! /usr/bin/env ruby
#
# Submit various scripts to PBS
#
# Somatic calling 
#
#   submit.rb somatic caller listfile
#
# Where caller is a script named $BASE/update_caller.sh The listfile should
# contain paired bam files, which can have either a full path, or without the
# path but living in /data/mapping/cancer or specified with --datapath.
