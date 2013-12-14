#! /usr/bin/env ruby
#
# Create a list of paired somatic files taking the base path as input. The
# assumption is that base file names are almost the same, with one different
# character R or N for normal and T for tumor.

dir=ARGV.shift

raise 'Need a directory as first parameter' if not File.directory?(dir)

$stderr.print "Fetching all BAM names in #{dir}\n"
bamlist = `find #{dir} -type f -name '*.bam'`.strip.split("\n").sort.uniq

# ---- Now pair them - dropping the singletons
list=[]
first = nil
bamlist.each do |fn|
  name = File.basename(fn)
  if first
    second = name
    # p [first,second] # we test a pair
    if first.size == second.size
      # count 'T'
      c1 = first.scan(/\w/).inject(Hash.new(0)){|h, c| h[c] += 1; h}
      c2 = second.scan(/\w/).inject(Hash.new(0)){|h, c| h[c] += 1; h}
      count_t = c2['T']-c1['T'] 
      if count_t == 1
        list << [first,second]
        first = nil
        next
      end
    end
  end
  first = name
end

list.each do |row|
  puts row.join("\t")
end


