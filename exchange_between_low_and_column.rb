#!/usr/bin/ruby

if ARGV.size == 0
  "print exchange_between_low_and_column.rb out.clst"
  exit
end

rows_h   = {}
columns_h = {}

a = []
ARGF.each_with_index do |x, i|
  a = x.chomp.split("\t")
  rows_h[ i ] = a
  a.each_with_index do |y, j|
    columns_h[j] = [] if columns_h[j] == nil
    columns_h[j][i] = y
  end
end

sorted_rows_keys    = rows_h.keys.sort
sorted_columns_keys = columns_h.keys.sort

sorted_columns_keys.each_with_index do |j, k|
  print columns_h[j].join("\t"), "\n"
end


