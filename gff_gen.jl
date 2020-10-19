function gen_gffs(gff::Vector{String},fa::String,vcf::String,win_exp::Int64,n_max::Int64)

    # Initialize VCF variables
    het_records = Vector{GFF3.Record}()
    hom_records = Vector{GFF3.Record}()
    reader_vcf = open(VCF.Reader,vcf)
    record = VCF.Record()
    read!(reader_vcf,record)

    # Initialize FASTA variables
    prev_end = 0
    reader_fa = open(FASTA.Reader, fa, index=fa*".fai")
    chr_names = reader_fa.index.names
    curr_chr = chr_names[1]
    fa_record = reader_fa[curr_chr]
    close(reader_fa)

    # Chrom sizes
    chr_sizes = reader_fa.index.lengths
    chr_size = chr_sizes[findfirst(x->x==curr_chr,chr_names)]

    # Loop over variants
    while true

        # Check if chromosome is in FASTA file. If not, obtain the next valid 1
        record = next_record(reader_vcf,chr_names,record)

        # Check we have loaded the right chromosome from FASTA
        if !(curr_chr == VCF.chrom(record))
            curr_chr = VCF.chrom(record)
            reader_fa = open(FASTA.Reader,fa,index=fa*".fai")
            fa_record = reader_fa[curr_chr]
            chr_size = chr_sizes[findfirst(x->x==curr_chr,chr_names)]
            close(reader_fa)
            prev_end=0
        end

        # Get entire (contiguous) haplotype using the PS field
        wSt = VCF.pos(record)
        wEnd = VCF.pos(record)
        het1 = Vector{Int64}()
        het2 = Vector{Int64}()
        if "PS" in VCF.format(record)
            # If PS tag, then phased
            ind_PS = findfirst(x->x=="PS",VCF.format(record))
            curr_ps = VCF.genotype(record)[1][ind_PS]
            record,wEnd = get_records_ps!(reader_vcf,fa_record,record,curr_ps,wEnd,het1,het2)
        else
            # If no PS tag, then single SNP (marked as NOPS)
            curr_ps = "NOPS"
            is_het_cpg!(record,fa_record,het1,het2)
            eof(reader_vcf) ? record=VCF.Record() : read!(reader_vcf, record)
        end
 # Expand window win_exp on each side
        win = expand_win(wSt,wEnd,win_exp,chr_size)

        # Obtain DNA sequence from reference and homozygous CpG loci from it
        wSeq = convert(String,FASTA.sequence(fa_record,win[1]:win[2]))
        cpg_pos = map(x->getfield(x,:offset),eachmatch(r"CG",wSeq)).+win[1].-1

        # Store heterozygous haplotype & corresponding CpG sites if homozygous CpG site/s
        if length(cpg_pos)>0

            # Obtain required number of models based on homozygous CpG sites
            n_mod = find_num_models(cpg_pos,n_max,1)
            mod_size = div(length(cpg_pos),n_mod)
            mod_rem = length(cpg_pos)-n_mod*mod_size

            # Print each set separately
            for i=1:n_mod
                # Index of subset of CpG sites in cpg_pos
                hom_rng = (mod_size*(i-1)+1):min((mod_size*i),length(cpg_pos))
                # Add remainders to last block
                i==n_mod && mod_rem>0 && (hom_rng=extrema(hom_rng)[1]:length(cpg_pos))
                # Get subset of CpG sites for each
                hom_mod = cpg_pos[hom_rng]
                # Get submodel boundaries
                inf_bound = i==1 ? win[1] : cpg_pos[max(1,extrema(hom_rng)[1]-1)]
                sup_bound = i==n_mod ? win[2] : hom_mod[end]
                # Get heterozygous CpG sites in current submodel
                het1_mod = het1[findall(x-> inf_bound <= x <= sup_bound,het1)]
                het2_mod = het2[findall(x-> inf_bound <= x <= sup_bound,het2)]
                # Check there is no overlap between homozygous and heterozygous CpG sites
                if !isempty(intersect(union(het1_mod,het2_mod),hom_mod))
                    print_log("Found overlap between homozygous and heterozygous CpG sites.")
                    print_log("Check FASTA file is masked. Exiting julia ...")
                    exit(1)
                end
                # Make submodel window a multiple of G (no new CpGs included)
                win_mod = expand_win(inf_bound,sup_bound,0,chr_size)
                # Push output string
                out_str = "$(curr_chr)\t.\t$(curr_ps)-$(i)/$(n_mod)\t$(win_mod[1])\t$(win_mod[2])"
                out_str *= "\t.\t.\t.\tN=$(length(hom_mod));CpGs=$(hom_mod)"
                length(het1_mod)>0 && (out_str*=";hetCpGg1=$(het1_mod)")
                length(het2_mod)>0 && (out_str*=";hetCpGg2=$(het2_mod)")
                # Push GFF record
             push!(het_records,GFF3.Record(out_str))
            end

        end

        # Check if het_records object is too big and empty it if so
        sizeof(het_records)>GFF_BUFFER && write_gff!(gff[1],het_records)

        # Obtain homozygous portion
        hom_win = [prev_end+1,win[1]-1]
        hom_win = [max(1,hom_win[1]),min(chr_size,hom_win[2])]

        # Obtain DNA sequence from reference and homozygous CpG loci from it
        wSeq = convert(String,FASTA.sequence(fa_record,hom_win[1]:hom_win[2]))
        cpg_pos = map(x->getfield(x,:offset),eachmatch(r"CG",wSeq)).+hom_win[1].-1

        # Store homozygous stretch of CpG sites found
        if length(cpg_pos)>0
            out_str = "$(curr_chr)\t.\t.\t$(hom_win[1])\t$(hom_win[2])\t$(length(cpg_pos))\t.\t.\t"
            out_str *= "CpGs=$(cpg_pos)"
            push!(hom_records,GFF3.Record(out_str))
        end

        # Update previous end
        prev_end = win[2]

        # Check if hom_records object is too big and empty it if so
        sizeof(hom_records)>GFF_BUFFER && write_gff!(gff[2],hom_records)

        # If new record empty break while loop
        VCF.haspos(record) || break

    end
    # Dump remaining records
    write_gff!(gff[1],het_records)
    write_gff!(gff[2],hom_records)

    # Return
    return nothing
