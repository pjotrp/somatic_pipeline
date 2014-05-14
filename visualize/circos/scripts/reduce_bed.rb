#! /usr/bin/env ruby
#
# Reduce multiple BED entries for the same exome/gene to one entry

h = {} ; hp1 = {} ; hp2 = {}
$stdin.each_line do | s |
  chr,p1,p2,name = s.split
  p1 = p1.to_i
  p2 = p2.to_i
  label = chr+'_'+name
  if not h[label] 
    h[label] = [chr,name]
    hp1[label] = p1
    hp2[label] = p2
  else
    hp1[label] = p1 if p1 < hp1[label]
    hp2[label] = p2 if p2 > hp2[label]
  end
end

h.each do | k,v |
  print [h[k][0],hp1[k],hp2[k],h[k][1]].join("\t"),"\n"
end
