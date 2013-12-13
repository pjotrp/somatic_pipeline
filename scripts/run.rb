#! /usr/bin/env ruby
#
# This runner takes a file of normal-tumor BAM files and executes the script
#
# Configuration happens in a JSON file
#
# E.g.
#
#   ./scripts/run.rb --config run.json ./scripts/varscan2.sh somatic_list.txt

require 'json'

def parse_args(args)
  options = { }

  consume = lambda { |args|
    return args if File.exist?(args[0]) # reached the executable command
    case args[0]
      when '--config', '-c'
        options[:config] = File.expand_path(args[1])
        consume.call(args[2..-1])
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
  # Write 'env.sh'
  File.open('env.sh','w') do | f |
    config.each do |k,v|
      print "config: ",k.to_s.upcase,"=\"",v,"\"\n"
      f.print k.to_s.upcase,"=\"",v,"\"\n"
    end
  end
end

File.read(listfn).each_line do | line |
  next if line =~ /^#/
  normal,tumor = line.strip.split(/\s+/)
  print "Normal=",normal,"\tTumor=",tumor
  normalname=File.basename(normal,'.bam')
  tumorname=File.basename(normal,'.bam')
  Kernel.system(["/bin/bash",script,(config ? '--config env.sh' : ''),normalname,tumorname,normal,tumor].join(" "))
end

