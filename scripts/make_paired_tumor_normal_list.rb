#! /usr/bin/env ruby
#
# Create a list of paired somatic files reading filenames as input. The
# assumption is that base file names are almost the same, with one different
# character R or N for normal and T for tumor. Here only T is checked.
#
# Example
#
#   ./make_paired_tumor_normal_list.rb /export/data/MBC/test/*p.bam
#   

# $stderr.print "Fetching all BAM names in #{dir}\n"
# bamlist = `find #{dir} -type f -name '*.bam'`.strip.split("\n").sort.uniq

bamlist = ARGV.sort.uniq

def match_pair(first,second)
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
        $stderr.print "WARNING: check match of #{first} with #{second}!\n"
        return nil
      end
      return [first,second]
    end
  end
  nil
end

# ---- Now pair them - dropping the singletons
list = []
first = nil
bamlist.each do |fn|
  if first
    second = fn
    res = match_pair(first,second)
    if res 
      list << res
      first = nil
      next
    end
  end
  $stderr.print "No match for #{first}\n" if first
  first = fn
end

list.each do |row|
  puts row.join(" ")
end



