import os

import geneplexus.cli


class PatchGP:
    avgps = None
    df_convert_out = None
    df_convert_out_subset = None
    df_probs = None
    df_edge = None
    df_edge_sym = None
    df_sim_GO = None
    df_sim_Dis = None

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
        skip_mdl_sim=True,
    )

    assert not os.path.isfile(geneplexus.cli.TMP_LOG_PATH)
