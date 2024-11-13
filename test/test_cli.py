import os

import geneplexus.cli


class PatchGP:
    avgps = None
    df_convert_out = None
    df_convert_out_subset = None
    df_probs = None
    df_edge = None
    df_edge_sym = None
    df_sim = None
    pos_genes_in_net = None
    negative_genes = None
    net_genes = None
    neutral_gene_info = None
    mdl_weights = None
    isolated_genes = None
    isolated_genes_sym = None

    def dump_config(self, path):
        return


def test_save_results(mocker, tmpdir):
    mocker.patch("numpy.savetxt")
    mocker.patch("geneplexus.cli.df_to_tsv")
    print(f"{tmpdir=}")
    geneplexus.cli.save_results(
        gp=PatchGP(),
        outdir=tmpdir,
        zip_output=False,
        overwrite=True,
        skip_mdl_sim=False,
        skip_sm_edgelist=False,
    )

    assert not os.path.isfile(geneplexus.cli.TMP_LOG_PATH)
