# Simple Somatic Pipeline

Various controller scripts for the SOLiD based somatic pipeline.

    ./scripts contains the pipeline control scripts

The scripts/update_xxxx are designed to run on their own, accepting
a name, normal and tumor BAM files, and perhaps some other arguments.

There is one control scripts which take a list of normal and tumor BAM
file pairs and a few other options. The control script can be used for
standalone running and for PBS.

For more information, see

    ./scripts/run.rb

One of the features of this pipeline is that is makes extensive use of 
[once-only](https://github.com/pjotrp/once-only).

## Copyright

Copyright (c) 2013,2014 Cuppen Group. See LICENSE.txt for further details.

