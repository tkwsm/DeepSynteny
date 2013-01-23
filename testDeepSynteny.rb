#!/usr/bin/ruby

require 'DeepSynteny.rb'
require 'runit/assert'
include RUNIT::Assert
include DeepSynteny

# Test for class MutualPair

  spe1gA = "36078"                  # scaffold_6
  spe1gA_position_from = 1508815
  spe1gA_position_to   = 1514256
  spe1gB = "215715"                 # scaffold_236
  spe1gB_position_from = 107644
  spe1gB_position_to   = 108317
  spe1gC = "228088"                 # scaffold_6
  spe1gC_position_from = 1516314
  spe1gC_position_to   = 1525322
  spe1gD = "238552"                 # scaffold_6
  spe1gD_position_from = 1555994
  spe1gD_position_to   = 1558281

  spe2gC = "aug_v2a.07489"          # scaf2807
  spe2gC_position_from = 38027
  spe2gC_position_to   = 50805
  spe2gD = "aug_v2a.19715"          # scaf11793
  spe2gD_position_from = 2848
  spe2gD_position_to   = 4715
  mp = MutualPair.new( spe1gC , spe2gC )
  assert_equal( spe1gC, mp.gene1 )
  assert_equal( spe2gC, mp.gene2 )

# Test for class ScaffoldDat

  species1_id        = "species1"
  spe1_scaf_id_1     = "scaffold_6"
  spe1_scaf_id_1_length = 2383264
  spe1_scaf_id_2     = "scaffold_236"
  spe1_scaf_id_2_length = 363062
  spe1sd1 = ScaffoldDat.new( species1_id , spe1_scaf_id_1, spe1_scaf_id_1_length )
  spe1sd2 = ScaffoldDat.new( species1_id , spe1_scaf_id_2, spe1_scaf_id_2_length )

  species2_id        = "species2"
  spe2_scaf_id_1     = "scaf2807"
  spe2_scaf_id_2     = "scaf11793"
  spe2_scaf_id_1_length = 58988
  spe2_scaf_id_2_length = 7649
  spe2sd1 = ScaffoldDat.new( species1_id , spe2_scaf_id_1, spe2_scaf_id_1_length )
  spe2sd2 = ScaffoldDat.new( species1_id , spe2_scaf_id_2, spe2_scaf_id_2_length )
 
  spe1sd1.add_gene( spe1gC )
  assert( spe1sd1.genes_of_the_scaffold.include?( spe1gC ), message="#{spe1gC} not included" )

# Test for class GeneDat

  spe1gCgd = GeneDat.new( spe1sd1, spe1gC_position_from, spe1gC_position_to, spe1sd1 )
  assert_equal( "scaffold_6", spe1gCgd.scaffold_id )

# Test for class GenesDat

  spe1gsd = GenesDat.new( species1_id )
  spe2gsd = GenesDat.new( species2_id )
  spe1gsd.add_new_gene( spe1gA, spe1gA_position_from, spe1gA_position_to, spe1sd1 )
  spe1gsd.add_new_gene( spe1gB, spe1gB_position_from, spe1gB_position_to, spe1sd2 )
  spe1gsd.add_new_gene( spe1gC, spe1gC_position_from, spe1gC_position_to, spe1sd1 )
  spe1gsd.add_new_gene( spe1gD, spe1gD_position_from, spe1gD_position_to, spe1sd1 )
  spe2gsd.add_new_gene( spe2gC, spe2gC_position_from, spe2gC_position_to, spe2sd1 )
  spe2gsd.add_new_gene( spe2gD, spe2gD_position_from, spe2gD_position_to, spe2sd1 )

# Test for class ScaffoldsDat

  spe1ssd = ScaffoldsDat.new( species1_id )
  spe1ssd.add_new_scaffold( species1_id, spe1_scaf_id_1, spe1_scaf_id_1_length)
  spe1ssd.add_new_scaffold( species1_id, spe1_scaf_id_2, spe1_scaf_id_2_length)
  spe2ssd = ScaffoldsDat.new( species2_id )

# Test for class MutualBestPoint

  mbp = MutualBestPoint.new( species1_id, species2_id, spe1ssd, spe2ssd, spe1gsd, spe2gsd )
  
  mbp.get_mutual_bests( spe1gC, spe2gC )
  mbp.get_mutual_bests( spe1gD, spe2gD )
  mutualbestid_C = mbp.spe1gene2mutualbestid( spe1gC )
  assert_equal( spe1gC, mbp.mutual_bests( mutualbestid_C ).gene1.geneid )

  spe1ssd.sort_by_scaffold_length!
  assert_equal( ["scaffold_6", "scaffold_236"], spe1ssd.scaffolds_order )
  spe1ssd.sort_reverse_by_scaffold_length!
  assert_equal( ["scaffold_236", "scaffold_6"], spe1ssd.scaffolds_order )

  print "\nAll Test for class DeepSynteny was cleared\n\n"
