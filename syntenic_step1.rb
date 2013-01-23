#!/usr/bin/ruby
require 'rubygems'
require 'bio'


if __FILE__ == $0

  if ARGV.size != 6
    print "\n USAGE: ruby syntenic_step1.rb < species1_gff > < species1_scaffolds > < species2_gff > < species2_scaffolds > < mutual_best_list >  < output_prefix > \n\n"
    exit
  else
  end

  require './DeepSynteny.rb'

  species1gff_file      = open( ARGV.shift )
  spe1scaffolds_file    = open( ARGV.shift )
  species2gff_file      = open( ARGV.shift )
  spe2scaffolds_file    = open( ARGV.shift )
  mutual_best_list_file = open( ARGV.shift )
  outprefix             = ARGV.shift
  output1      = File.new("#{outprefix}.plot", "w+")
  output2      = File.new("#{outprefix}.x", "w+")
  output3      = File.new("#{outprefix}.y", "w+")
  output4      = File.new("#{outprefix}.clst", "w+")

  species1_name = "species_1"
  species2_name = "species_2"

  mutualbests = DeepSynteny::easy_pipeline_1( species1_name, species2_name, spe1scaffolds_file, spe2scaffolds_file, species1gff_file, species2gff_file, mutual_best_list_file )

# OUTPUT mutual bests plot (PREFIX.plot)

  mutualbests.each do |mutual_best|
    output1.print mutual_best.gene1.geneid, "\t"
    output1.print mutual_best.gene2.geneid, "\t"
    output1.print mutual_best.gene1.coordinate_position_from, "\t"
    output1.print mutual_best.gene2.coordinate_position_from, "\n"
  end

# OUTPUT x axis lines of scaffolds (PREFIX.x)

  basic_data = 0
  mutualbests.spe1scaf_obj.each_scaffold do |sobj|
    basic_data += sobj.scaffold_length
    output2.print sobj.scaffold_id, "\t"
    output2.print sobj.genes_of_the_scaffold.size, "\t"
    output2.print sobj.scaffold_length, "\t"
    output2.print basic_data, "\n"
  end

# OUTPUT y axis lines of scaffolds (PREFIX.y)
  basic_data = 0
  mutualbests.spe2scaf_obj.each_scaffold do |sobj|
    basic_data += sobj.scaffold_length
    output3.print sobj.scaffold_id, "\t"
    output3.print sobj.genes_of_the_scaffold.size, "\t"
    output3.print sobj.scaffold_length, "\t"
    output3.print basic_data, "\n"
  end

# OUTPUT mutual best matrix fo cluster (PREFIX.clust)

  spe1s_to_spe2s_h = {}
  spe2s_to_spe1s_h = {}
  mutualbestid = ""
  spe1sca = ""
  spe2sca = ""
  mutualbests.spe1scaf_obj.each_scaffold do |sobj1|
    sobj1.genes_of_the_scaffold.each do |gene1|
       if mutualbests.spe1gene2mutualbestid( gene1 )
         mutualbestid = mutualbests.spe1gene2mutualbestid( gene1 )
         spe1s = mutualbests.mutual_bests( mutualbestid ).gene1.scaffold_id
         spe2s = mutualbests.mutual_bests( mutualbestid ).gene2.scaffold_id

         spe1s_to_spe2s_h[ spe1s ] = {} if spe1s_to_spe2s_h[ spe1s ] == nil
         spe1s_to_spe2s_h[ spe1s ][ spe2s ] = 0 if spe1s_to_spe2s_h[ spe1s ][ spe2s ] == nil
         spe1s_to_spe2s_h[ spe1s ][ spe2s ] += 1

         spe2s_to_spe1s_h[ spe2s ] = {} if spe2s_to_spe1s_h[ spe2s ] == nil
         spe2s_to_spe1s_h[ spe2s ][ spe1s ] = 0 if spe2s_to_spe1s_h[ spe2s ][ spe1s ] == nil
         spe2s_to_spe1s_h[ spe2s ][ spe1s ] += 1
       end
    end
  end

  spe1s_keys = spe1s_to_spe2s_h.keys
  spe2s_keys = spe2s_to_spe1s_h.keys

  output4.print "spe2-spe1"
  spe1s_keys.each{|spe1s| output4.print "\t#{spe1s}" }
  output4.print "\n"

  spe2s_keys.each do |spe2s|
    output4.print "#{spe2s}"
    spe1s_keys.each do |spe1s|
      val = 0.0
      if spe2s_to_spe1s_h[ spe2s ][ spe1s ]
        val = spe2s_to_spe1s_h[ spe2s ][ spe1s ].to_i
      else
      end
      output4.print "\t#{val}"
    end
    output4.print "\n"
  end

##########################################################

end


