#! /usr/bin/env ruby

require 'csv'
require 'solid_assert'

FREQ=8

SolidAssert.enable_assertions
h=ARGV.map{ |s| s.split('=') }.to_h
p h

csv_parse = lambda { |cmd,args = "NONE=0"|
  fcmd = "env #{args} #{ENV['HOME']}/izip/git/opensource/ruby/bioruby-rdf/scripts/sparql-csv.sh "+cmd
  $stderr.print "===> Calling: #{fcmd}\n"
  CSV::parse(`#{fcmd}`).drop(1)
}

no_parse = lambda { |cmd,args = "NONE=0"|
  fcmd = "env #{args} #{ENV['HOME']}/izip/git/opensource/ruby/bioruby-rdf/scripts/sparql-csv.sh "+cmd
  $stderr.print "===> Calling: #{fcmd}\n"
  `#{fcmd}`
}


samples = csv_parse.call("../count_samples.rq").flatten

p samples

Dir.mkdir("output") if not File.exist?("output")

count = {}
ref   = {}
samples.each do | s |
  File.open("output/"+s+'.tsv','w') do | f |
    res = no_parse.call("germline.rq","SAMPLE=#{s}")
    table = CSV::parse(res).drop(1)
    t_has = {}
    f.print "gene_name,chr,pos,alt,sample,is_cancer,is_bc,is_ova,freq,dbsnp,cosmic,remark,info,type".split(/,/).join("\t"),"\n"
    table.each do | row |
      obj = row[0..3]
      cfreq = row[FREQ].to_f
      if not t_has[obj]
        t_has[obj] = true
        count[obj] = 0 if !count[obj]
        count[obj] += 1
        ref[obj] = row
      else
        # plug the largest frequency
        freq = ref[obj][FREQ].to_f
        ref[obj][FREQ] = cfreq if freq == nil or freq<cfreq
      end
    end
    ref.each do | k,row |
      freq = row[FREQ].to_f
      if freq and freq <= 0.10 
        f.print row.join("\t"),"\n"
      end
    end
    # f.print res
  end
end

File.open("output/germline.tsv",'w') do | f |
  f.print "num,gene,chr,pos,alt,sample,is_cancer,breast,ovarian,freq,dbsnp,cosmic,remark,info,type".split(/,/).join("\t"),"\n"
  count.sort_by{|k,v| v}.reverse.each { | k,v | 
    freq = ref[k][FREQ].to_f
    if freq and freq <= 0.10
      f.print v,"\t",ref[k].join("\t"),"\n" if v > 1
    end
  }
end
