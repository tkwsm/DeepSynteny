#!/usr/bin/ruby
require 'rubygems'
require 'bio'

if __FILE__ == $0

  if ARGV.size != 8
    print "\n ruby syntenic_step2.rb < species1_gff > < species1_scaffolds > < species2_gff > < species2_scaffolds > < mutual_best_list > < species1_hint_file > < species2_hint_file > <output_prefix> \n\n"
    exit
  else
  end

#  require '/genefs/SatohU/takeshik/synteny/n_vectensis_vs_a_digitifera/mututal_best_dat.rb'
  require './DeepSynteny.rb'

  species1gff_file      = open( ARGV.shift )
  spe1scaffolds_file    = open( ARGV.shift )
  species2gff_file      = open( ARGV.shift )
  spe2scaffolds_file    = open( ARGV.shift )
  mutual_best_list_file = open( ARGV.shift )
  species1_hint_file    = open( ARGV.shift )
  species2_hint_file    = open( ARGV.shift )
  outprefix             = ARGV.shift
  output1      = File.new("#{outprefix}.2.plot", "w+")
  output2      = File.new("#{outprefix}.2.x", "w+")
  output3      = File.new("#{outprefix}.2.y", "w+")

  species1_name = "species_1"
  species2_name = "species_2"

  mutualbests = DeepSynteny::easy_pipeline_2( species1_name, species2_name, spe1scaffolds_file, spe2scaffolds_file, species1gff_file, species2gff_file, mutual_best_list_file, species1_hint_file, species2_hint_file )

  mutualbests.each do |mutual_best|
    output1.print mutual_best.gene1.geneid, "\t"
    output1.print mutual_best.gene2.geneid, "\t"
    output1.print mutual_best.gene1.coordinate_position_from, "\t"
    output1.print mutual_best.gene2.coordinate_position_from, "\n"
  end

  basic_data = 0
  mutualbests.spe1scaf_obj.each_scaffold do |sobj|
    basic_data += sobj.scaffold_length
    output2.print sobj.scaffold_id, "\t"
    output2.print sobj.genes_of_the_scaffold.size, "\t"
    output2.print sobj.scaffold_length, "\t"
    output2.print basic_data, "\n"
  end

  basic_data = 0
  mutualbests.spe2scaf_obj.each_scaffold do |sobj|
    basic_data += sobj.scaffold_length
    output3.print sobj.scaffold_id, "\t"
    output3.print sobj.genes_of_the_scaffold.size, "\t"
    output3.print sobj.scaffold_length, "\t"
    output3.print basic_data, "\n"
  end

end


