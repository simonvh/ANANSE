import ananse.exceptions
import ananse.utils
import ananse.enhancer


def enhancer(**kwargs):
    """
    CLI parser for ananse.enhancer.Enhancer
    """
    # print(dict(**kwargs))
    # print(kwargs.get("genome"))
    genome = kwargs.get("genome", None)
    genome = ananse.utils.check_genome(genome)

    # etype = kwargs.get("type")
    # etype = ananse.utils.check_type(etype)
    #
    # bam = kwargs.get("bam")
    # bam = ananse.utils.check_bam(bam)

    peak = kwargs.get("peak")
    peak = ananse.utils.check_file(peak, "peak")

    output = kwargs.get("output")
    output = ananse.utils.check_output(output, "enhancer.bed")

    signal_column = kwargs.get("signal")
    summit_column = kwargs.get("summit")

    b = ananse.enhancer.Enhancer(
        genome=genome,
        #bam=bam,
        peak=peak,
        output=output,
        signal_column=signal_column,
        summit_column=summit_column
    )
    b.run_enhancer()
        #bam, peak, output
    #)

    # # if not os.path.exists(args.fin_rpkm):
    # #     print("File %s does not exist!" % args.fin_rpkm)
    # #     sys.exit(1)
    # genome=args.genome
    # etype=args.etype
    #
    # if genome == "hg38" and etype == "hg38H3K27ac":
    #     b = ananse.enhancer.Enhancer(
    #         genome=args.genome,
    #         bam=args.bam,
    #         epeak=args.epeak,
    #         bed_output=args.bed_output
    #     )
    #     b.run_enhancer(
    #         args.bam, args.epeak, args.bed_output
    #     )
    # elif etype == "p300":
    #     b = ananse.enhancer.P300Enhancer(
    #         genome=args.genome,
    #         bam=args.bam,
    #         epeak=args.epeak,
    #         bed_output=args.bed_output
    #     )
    #     b.run_enhancer(
    #         args.bam, args.epeak, args.bed_output
    #     )
    # elif etype == "ATAC":
    #     b = ananse.enhancer.AtacEnhancer(
    #         genome=args.genome,
    #         bam=args.bam,
    #         epeak=args.epeak,
    #         bed_output=args.bed_output
    #     )
    #     b.run_enhancer(
    #         args.bam, args.epeak, args.bed_output
    #     )
