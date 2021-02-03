import ananse.utils
import ananse.binding


def binding(**kwargs):
    """
    CLI parser for ananse.binding.Binding
    """
    genome = kwargs.get("genome", None)
    genome = ananse.utils.check_genome(genome)

    enhancer = kwargs.get("enhancers", None)
    enhancer = ananse.utils.check_file(enhancer, "enhancers")

    # TODO: check for both required files
    pfm = kwargs.get("pfm")
    #pfm = ananse.utils.check_pfm(pfm)

    ncpu = kwargs.get("ncpu")

    output = kwargs.get("output")
    output = ananse.utils.check_output(output, "binding.bed")

    a = ananse.binding.Binding(
        ncore=ncpu,
        genome=genome,
        # gene_bed=args.annotation,
        pfmfile=pfm,
        #include_notfs=args.include_notfs,
        #rm_curated=args.rm_curated,
        #etype=args.etype,
        #tffile=args.tffile
    )
    a.run_binding(enhancer, output)

# #!/usr/bin/env python
# # Copyright (c) 2009-2019 Quan Xu <qxuchn@gmail.com>
# #
# # This module is free software. You can redistribute it and/or modify it under
# # the terms of the MIT License, see the file COPYING included with this
# # distribution.
#
# from __future__ import print_function
# import sys
# import os
#
# import ananse.binding
#
#
# def binding(args):
#     if not os.path.exists(args.fin_rpkm):
#         print("File %s does not exist!" % args.fin_rpkm)
#         sys.exit(1)
#
#     a = ananse.binding.Binding(
#         ncore=args.ncore,
#         genome=args.genome,
#         # gene_bed=args.annotation,
#         pfmfile=args.pfmfile,
#         include_notfs=args.include_notfs,
#         rm_curated=args.rm_curated,
#         etype=args.etype,
#         tffile=args.tffile
#     )
#     a.run_binding(args.fin_rpkm, args.outfile)
