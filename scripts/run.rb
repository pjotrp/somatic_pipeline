#! /usr/bin/env ruby
#
# This runner takes a file of normal-tumor BAM files and executes the script standalone
# or on PBS with the --pbs switch.
#
# Configuration happens in a JSON file
#
# E.g.
#
#   ./scripts/run.rb --config run.json ./scripts/varscan2.sh somatic_list.txt
#
# or a full on with --pbs
#
#   wgs01:~/data/trials/chr17$ ~/opt/somatic_pipeline/scripts/run.rb --pbs --config run.json ~/opt/somatic_pipeline/scripts/varscan2.sh somatic2_bams.txt
#
# Example of run.json
#
# {
#   "refseq" : "/export/data/MBC/test/GRCh37_gatk.fasta",
#   "samtools": "/home/wrk/opt/bin/samtools",
#   "dataroot": "/export/data/MBC"
# }
#   

require 'fileutils'
require 'json'

if ARGV.size == 0
  print USAGE
  exit 1
end

def exit_error errval = 1, msg = nil
  $stderr.print msg if msg
  $stderr.print "\n**ERROR** once-only returned error #{errval}\n"
  exit errval
end

def parse_args(args)
  options = { }

  consume = lambda { |args|
    return args if File.exist?(args[0]) # reached the executable command
    case args[0]
      when '--config', '-c'
        options[:config] = File.expand_path(args[1])
        consume.call(args[2..-1])
      when '--pbs'
        options[:pbs] = args[1]
        consume.call(args[1..-1])
      when '--first','-1'
        options[:first] = true
        consume.call(args[1..-1])
      else
        $stderr.print "**ERROR** Can not parse arguments",args
        exit_error(1)
      end
  }

  return consume.call(args),options
end

args,options = parse_args(ARGV)

script = args[0]
raise "Expected a valid script" if not File.executable?(script)
listfn = args[1]
raise "Expected a list file" if not listfn

config = if options[:config]
           json = File.read(options[:config])
           JSON.parse(json,:symbolize_names => true)
         end

if config
  # ---- Use the JSON options and write 'env.sh'
  env_sh = 'env_'+File.basename(script)
  File.open(env_sh,'w') do | f |
    config.each do |k,v|
      print "config: ",k.to_s.upcase,"=\"",v,"\"\n"
      f.print k.to_s.upcase,"=\"",v,"\"\n"
    end
  end
end

abs_env_sh = File.absolute_path(env_sh)

File.read(listfn).each_line do | line |
  next if line =~ /^#/
  normal,tumor = line.strip.split(/\s+/)
  find = lambda { |bam|
    # ---- Find the full path of the file name in the list
    res ||= `find #{config[:dataroot]} -name #{bam}`.strip
    raise "Too many candidates for #{bam}" if res.split(/\n/).size > 1
    # ---- Make the file read-only if we can write to it
    File.chmod(0444,res)
    # ---- From now on use a symlink to the file
    FileUtils.ln_s(res,bam,verbose: true)
    bam
  }
  normal = find.call(normal) if not File.exist?(normal)
  tumor = find.call(tumor) if not File.exist?(tumor)
  print "Normal=",normal,"\tTumor=",tumor
  normalname=File.basename(normal,'.bam')
  tumorname=File.basename(tumor,'.bam')
  p [normalname,tumorname]
  cmd = [script,(config ? '--config '+abs_env_sh : ''),normalname,tumorname,normal,tumor].join(" ")
  if options[:pbs]
    # ---- Submit to PBS
    p cmd
    jobname=tumorname
    jobname='mbc'+tumorname if jobname =~ /^\d/
    print `echo \"#{cmd}\" | qsub -P SAP42 -N #{jobname} -cwd`
  else
    # ---- Run standalone
    p cmd
    Kernel.system("/bin/bash "+cmd)
    if $?.exitstatus 
      $stderr.print "Command <"+cmd+"> did not complete!"
      exit_error(1)
    end
  end
  if options[:first]
    $stderr.print "Stopped after --first"
    exit_error(1)
  end
end

