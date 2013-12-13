Various controllers for the somatic pipeline

./runners contains little programs used in a pipeline
./scripts contains the pipeline control scripts

The scripts/update_xxxx are designed to run on their own, accepting
a name, normal and tumor BAM files, and perhaps some other arguments.

There are two control scripts which take a list of normal and tumor
BAM file pairs and a few other options. One control script is for
standalone running (run.rb) and the other is for PBS (submit.rb).

For more information type

./scripts/run.rb

or

./scripts/submit.rb

