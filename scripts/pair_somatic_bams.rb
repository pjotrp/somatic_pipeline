#! /usr/bin/env ruby
#
# Create a list of paired somatic files taking the base path as input. The
# assumption is that base file names are almost the same, with one different
# character R or N for normal and T for tumor. Here only T is checked.

dir=ARGV.shift

raise 'Need a directory as first parameter' if not File.directory?(dir)

$stderr.print "Fetching all BAM names in #{dir}\n"
bamlist = `find #{dir} -type f -name '*.bam'`.strip.split("\n").sort.uniq

# ---- Now pair them - dropping the singletons
list=[]
first = nil
bamlist.each do |fn|
  if first
    second = fn
    # p [first,second] # we test a pair
    if first.size == second.size
      b1 = File.basename(first)
      b2 = File.basename(second)
      # count 'T'
      c1 = b1.scan(/\w/).inject(Hash.new(0)){|h, c| h[c] += 1; h}
      c2 = b2.scan(/\w/).inject(Hash.new(0)){|h, c| h[c] += 1; h}
      count_t = c2['T']-c1['T'] 
      if count_t > 0
        # Final comparison
        if b1.sub(/R/,'T') != b2 and b1.sub(/N/,'T') != b2
          $stderr.print "WARNING: check match of #{b1} with #{b2}!\n"
        end
        list << [first,second]
        first = nil
        next
      end
    end
  end
  first = fn
end

list.each do |row|
  puts row.join(" ")
end



