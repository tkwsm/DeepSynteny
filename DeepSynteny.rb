#!/usr/bin/ruby
require 'rubygems'
require 'bio'

module DeepSynteny

  class MutualPair

    def initialize( species1geneA, species2geneA)
      @gene1 = species1geneA
      @gene2 = species2geneA
    end

    attr_reader :gene1, :gene2

  end

  class ScaffoldDat

    def initialize( species_id, scaffold_id, scaffold_length )
      if scaffold_length.class == String and scaffold_length !~ /\d/
        print "\nError! Scaffold_length should be integer!\n\n"; exit ;
      end
      @species_id      = species_id
      @scaffold_id     = scaffold_id
      @scaffold_length = scaffold_length.to_i
      @sorting_order   = 1
      @genes_of_the_scaffold = []
    end

    attr_reader   :scaffold_id, :scaffold_length, :species_id
    attr_accessor :sorting_order, :genes_of_the_scaffold

    alias :genes :genes_of_the_scaffold

    def add_gene( gene_id )
      @genes_of_the_scaffold << gene_id
    end

    def sort_gene_by_location!( genes_obj )
      @genes_of_the_scaffold.sort!{|x, y| genes_obj.info(x).position_from <=> genes_obj.info(y).position_from }
    end

  end

  class ScaffoldsDat

    def initialize( species_id )
      @species_id = species_id
      @scaffolds  = {}
      @scaffolds_with_mutual_best  = {}
      @scaffolds_order = []
    end

    attr_reader :scaffolds_order, :scaffolds_with_mutual_best

    def recreate_scaffolds!

      tmp_keys = @scaffolds.keys
      tmp_keys.each do |scaffold_id|

        if @scaffolds_with_mutual_best[ scaffold_id ] == nil
          @scaffolds.delete( scaffold_id ) 
          @scaffolds_order.delete( scaffold_id )
        elsif @scaffolds_with_mutual_best[ scaffold_id ].size == 0
          @scaffolds.delete( scaffold_id ) 
          @scaffolds_order.delete( scaffold_id )
        end

      end

    end

    def scaffolds_with_mutual_best( scaffold_id )
      @scaffolds_with_mutual_best[ scaffold_id ]
    end

    def add_scaffolds_with_mutual_best( sid, mutualbest_id ) #sid == scaffold_id
      @scaffolds_with_mutual_best[ sid ] = [] if @scaffolds_with_mutual_best[ sid ] == nil
      @scaffolds_with_mutual_best[ sid ] << mutualbest_id
      @scaffolds_order << sid unless @scaffolds_order.include?( sid )
    end

    def add_new_scaffold( species_id, scaffold_id, scaffold_length )
      scafdat = ScaffoldDat.new( species_id, scaffold_id, scaffold_length )
      @scaffolds[ scaffold_id ] = scafdat
      @scaffolds_order << scaffold_id unless @scaffolds_order.include?( scaffold_id )
    end

    def add_new_scaffold_longer_than( species_id, scaffold_id, scaffold_length, threshold_length )
      if scaffold_length > threshold_length
         add_new_scaffold( species_id, scaffold_id, scaffold_length )
      end
    end

    def info( scaffold_id )
      @scaffolds[ scaffold_id ]
    end

    def each_key
      for x in @scaffolds.keys
        yield x
      end
    end

    def each_scaffold
      @scaffolds_order.each do |scaffold_id|
        yield @scaffolds[ scaffold_id ]
      end
    end

    def each
      for k in @scaffolds.keys
        yield @scaffolds[ k ]
      end
    end

    def sort_by_total_genes!
      @scaffolds_order = @scaffolds.keys.sort!{|x, y| self.info(y).genes.size <=> self.info(x).genes.size }
      @scaffolds_order
    end

    def sort_by_scaffold_length!
      @scaffolds_order = @scaffolds_order.sort!{|x, y| info(y).scaffold_length <=> info(x).scaffold_length }
      @scaffolds_order
    end

    def sort_reverse_by_scaffold_length!
      @scaffolds_order = @scaffolds_order.sort!{|x, y| info(x).scaffold_length <=> info(y).scaffold_length }
      @scaffolds_order
    end

    def sort_by_hint!( hint_hash )
      @scaffolds.keys.each do |k|
        hint_hash[ k ] = 0 if hint_hash[ k ] == nil
      end
      @scaffolds_order = @scaffolds.keys.sort!{|x, y| hint_hash[y].to_i <=> hint_hash[x].to_i }
      @scaffolds_order
    end

    def sort_gene_by_location!( genes_obj )
      @scaffolds.each_key do |k|
        @scaffolds[k].sort_gene_by_location!( genes_obj )
      end
    end

  end

  class GeneDat

    def initialize( geneid, position_from, position_to, scaffold_obj )
      @geneid                   = geneid
      @scaffold_id              = scaffold_obj.scaffold_id
      @scaffold_length          = scaffold_obj.scaffold_length
      @position_from            = position_from
      @position_to              = position_to
      @coordinate_position_from = @position_from
      @coordinate_position_to   = @position_to
      scaffold_obj.add_gene( geneid )
    end

    attr_accessor :scaffold_id, :scaffold_length, :geneid,
                  :position_from, 
                  :position_to,
                  :coordinate_position_from, 
                  :coordinate_position_to

    def coord_from=( pos )
      @coordinate_position_from = pos
    end

  end

  class GenesDat

    def initialize( species_id )
      @species_id = species_id
      @genes      = {}
    end

    def add_new_gene( geneid, position_from, position_to, scaffold_obj )
      genedat = GeneDat.new( geneid, position_from, position_to, scaffold_obj )
      @genes[ geneid ] = genedat
    end

    def info( geneid )
      @genes[ geneid ]
    end

    def each_key
      for x in @genes.keys
        yield x
      end
    end

  end

  class MutualBestPoint

    def initialize( species1, species2, species1scaffolds_obj, species2scaffolds_obj, species1genes_obj, species2genes_obj )
      @mutual_bests   = {}
      @species1 = species1
      @species2 = species2
      @spe1scaf_obj  = species1scaffolds_obj
      @spe2scaf_obj  = species2scaffolds_obj
      @spe1genes_obj = species1genes_obj
      @spe2genes_obj = species2genes_obj
      @spe1_gene_coordinate = DeepSynteny::get_genes_coordinate( species1, @spe1scaf_obj, @spe1genes_obj )
      @spe2_gene_coordinate = DeepSynteny::get_genes_coordinate( species2, @spe2scaf_obj, @spe2genes_obj )
      @spe1gene2mutualbestids = {}
      @spe2gene2mutualbestids = {}
    end

    attr_accessor :spe1scaf_obj, :spe2scaf_obj

    attr_reader   :spe1gene2mutualbestids, :spe2gene2mutualbestids

    def mutual_bests( mutualbestid )
      @mutual_bests[ mutualbestid ]
    end

    def spe1gene2mutualbestid( gid1 )
      @spe1gene2mutualbestids[ gid1 ]
    end

    def get_mutual_bests( gid1, gid2)
      if @spe1genes_obj.info( gid1 ) and @spe2genes_obj.info( gid2 )
        mutualbestid = @mutual_bests.keys.size + 1
        @mutual_bests[ mutualbestid ] = MutualPair.new( @spe1genes_obj.info( gid1 ), @spe2genes_obj.info( gid2 ) )
        @spe1gene2mutualbestids[ gid1 ] = mutualbestid
        @spe2gene2mutualbestids[ gid2 ] = mutualbestid
        spe1_scaffoldid = @spe1genes_obj.info( gid1 ).scaffold_id
        spe2_scaffoldid = @spe2genes_obj.info( gid2 ).scaffold_id
        @spe1scaf_obj.add_scaffolds_with_mutual_best( spe1_scaffoldid, mutualbestid) 
        @spe2scaf_obj.add_scaffolds_with_mutual_best( spe2_scaffoldid, mutualbestid) 
      end
    end

    def get_mutual_bests_from_file( mutual_best_list_file )
      a = []
      mutual_best_list_file.each do |x|
        a = x.chomp.split("\t")
        next if x !~ /\S+/
        next if a.size < 2
        gid1 = a[0] 
        gid2 = a[1] 
        get_mutual_bests( gid1, gid2 )
      end
    end

    def sort_by_hint!( hint_hash, species )
      if    species == @species1
        @spe1scaf_obj.sort_by_hint!( hint_hash )
      elsif species == @species2
        @spe2scaf_obj.sort_by_hint!( hint_hash )
      else
        STDER.print "species should be #{@species1} or #{@species2}\n"
        exit
      end
      @spe1_gene_coordinate = DeepSynteny::get_genes_coordinate( @species1, @spe1scaf_obj, @spe1genes_obj )
      @spe2_gene_coordinate = DeepSynteny::get_genes_coordinate( @species2, @spe2scaf_obj, @spe2genes_obj )
    end

    def recreate_gene_coodinate!
      @spe1_gene_coordinate = DeepSynteny::get_genes_coordinate( @species1, @spe1scaf_obj, @spe1genes_obj )
      @spe2_gene_coordinate = DeepSynteny::get_genes_coordinate( @species2, @spe2scaf_obj, @spe2genes_obj )
    end

    def each_key
      for k in @mutual_bests.keys
        yield k
      end
    end

    def each
      for k in @mutual_bests.keys
        yield @mutual_bests[ k ]
      end
    end

    def each_pair
      for k in @mutual_bests.keys
        yield @mutual_bests[ k ]
      end
    end

  end

  def DeepSynteny::get_scaffolds_object( sca_fasta, species_id, min_scaf )
    minimum_length_for_scaffold = min_scaf
    scaffolds_object = ScaffoldsDat.new( species_id )
    scaff = Bio::FlatFile.new(Bio::FastaFormat, sca_fasta )
    scaff.each do |e|
      scaffold_id = e.definition
      scaffold_length = e.naseq.size
#      scaffolds_object.add_new_scaffold( species_id, scaffold_id, scaffold_length )
      scaffolds_object.add_new_scaffold_longer_than( species_id, scaffold_id, scaffold_length, minimum_length_for_scaffold )
    end
    return scaffolds_object
  end

  def DeepSynteny::get_genes_object( gff, species_id, scaffolds_object )

    genesdat     = GenesDat.new( species_id )

    a = []
    type = ""
    species_id = species_id
    scaffold_id = ""
    scaffold_length = 0
    gid = ""

    gff.each do |x|
      a = x.chomp.split("\t")
      type = a[2]
      next if type != "gene"
      if a[8] =~ /PID/
        gid = a[8].slice(/PID=([^\;]+)/, 1) #
      else
        gid = a[8].slice(/ID=([^\;]+)/, 1) #
      end
      pos_from = a[3].to_i
      pos_to   = a[4].to_i
      scaffold_id     = a[0].slice(/^(\S+)/, 1)            #
      scaffold_obj = scaffolds_object.info( scaffold_id )
      next if scaffold_obj == nil
      genesdat.add_new_gene( gid, pos_from, pos_to, scaffold_obj )
    end
    return genesdat
  end

  def DeepSynteny::get_genes_coordinate( species_id, scaffolds_obj, genes_obj  )

    current_x_axis_base = 0
    scaffolds_obj.sort_gene_by_location!( genes_obj )
    scaffolds_obj.scaffolds_order.each_with_index do |scaffold_id, i|
      next if scaffolds_obj.info( scaffold_id ) == nil

      scaffold_length = scaffolds_obj.info( scaffold_id ).scaffold_length
      scaffolds_obj.info( scaffold_id ).genes.each do |geneid|
        position_from = genes_obj.info(geneid).position_from.to_i
        genes_obj.info(geneid).coord_from = current_x_axis_base + position_from
      end
      current_x_axis_base += scaffold_length.to_i
    end
  end

  def DeepSynteny::get_cutree_file( cutree_output )
    output_h = {}
    a = []
    scafid = ""
    scaf_cluster_id = ""
    cutree_output.each do |x|
      a = x.chomp.split("\s")
      if scafid =~ /^\"/
        scafid = a[0].slice(/\"(\S+)\"/, 1)
      else
        scafid = a[0].slice(/(\S+)/, 1)
      end
      scaf_cluster_id = a[1].to_i
      output_h[ scafid ] = scaf_cluster_id
    end
    return output_h
  end

  def DeepSynteny::get_scaffold_cluster( scaffold_order_file )
    a = []
    s = ""
    c = ""
    h = {}
    scaffold_order_file.each_with_index do |x, i|
      a = x.chomp.split("\s")
      s += a.join("\t") if ( i % 2 == 0)
      c += a.join("\t") if ( i % 2 == 1)
    end
    c.chomp.split("\t").each_with_index do |cn, j|
      sn = s.chomp.split("\t")[j]
      h[ sn ] = cn
    end
    return h
  end

  def DeepSynteny::easy_pipeline_1( species1_name, 
                                      species2_name,
                                      spe1scaffolds_file, 
                                      spe2scaffolds_file, 
                                      species1gff_file,
                                      species2gff_file, 
                                      mutual_best_list_file )

    min_scaf_for_species1 = 20000
    min_scaf_for_species2 = 20000
    spe1scaf_obj = DeepSynteny::get_scaffolds_object( spe1scaffolds_file, species1_name, min_scaf_for_species1 )
    spe2scaf_obj = DeepSynteny::get_scaffolds_object( spe2scaffolds_file, species2_name, min_scaf_for_species2 )
    spe1genes_obj = DeepSynteny::get_genes_object( species1gff_file, species1_name, spe1scaf_obj )
    spe2genes_obj = DeepSynteny::get_genes_object( species2gff_file, species2_name, spe2scaf_obj )
    mutualbests = DeepSynteny::MutualBestPoint.new( species1_name, species2_name, spe1scaf_obj, spe2scaf_obj, spe1genes_obj, spe2genes_obj )
    mutualbests.get_mutual_bests_from_file( mutual_best_list_file )
    mutualbests.spe1scaf_obj.recreate_scaffolds!
    mutualbests.spe2scaf_obj.recreate_scaffolds!
#    mutualbests.spe1scaf_obj.sort_by_total_genes!
#    mutualbests.spe2scaf_obj.sort_by_total_genes!
    mutualbests.spe1scaf_obj.sort_by_scaffold_length!
    mutualbests.spe2scaf_obj.sort_by_scaffold_length!
    mutualbests.spe1scaf_obj.sort_by_scaffold_length!
    mutualbests.spe2scaf_obj.sort_by_scaffold_length!
    mutualbests.recreate_gene_coodinate!

    return mutualbests
 
  end

  def DeepSynteny::easy_pipeline_2( species1_name, 
                                      species2_name,
                                      spe1scaffolds_file, 
                                      spe2scaffolds_file, 
                                      species1gff_file,
                                      species2gff_file, 
                                      mutual_best_list_file, 
                                      spe1hint_file,
                                      spe2hint_file )

    mutualbests = DeepSynteny::easy_pipeline_1( species1_name, 
                                      species2_name,
                                      spe1scaffolds_file, 
                                      spe2scaffolds_file, 
                                      species1gff_file,
                                      species2gff_file, 
                                      mutual_best_list_file )

    mutualbests.spe1scaf_obj.recreate_scaffolds!
    mutualbests.spe2scaf_obj.recreate_scaffolds!
    mutualbests.spe1scaf_obj.sort_by_hint!( spe1hint )
    mutualbests.spe2scaf_obj.sort_by_hint!( spe2hint )
    mutualbests.recreate_gene_coodinate!

    return mutualbests
 
  end

end

