HOUR = 60
DAY = 24 * HOUR

parameters = {
    "vinf_incr":   0.25  / HOUR,
    "vrna_incr":   0.5   / HOUR,
    "vprot_incr":  0.167 / HOUR,

    "pirf3_incr":  0.75  / HOUR,
    "pirf3_decr":  0.125 / HOUR,

    "ifni_incr":   0.75  / HOUR,
    "ifni_decr":   0.25  / HOUR,

    "mm_pstat":   500.0,
    "pstat_incr":  40.0  / HOUR,
    "pstat_decr":  10.0  / HOUR,

    "isg_incr":    0.10  / HOUR,
    "isg_decr":    0.033 / HOUR,

    "q_ifne":      1.0   /  DAY,
    "k_ifn_sec":   5e5   /  HOUR,

    "vprot_inh_ifni":    2.0,
    "vprot_inh_pirf3":   2.0,
    "vprot_inh_pstat":   1.5,

    "isg_inh_vrna":      2.0,
    "isg_inh_vprot":     2.0,

    # zeroed transitions
    "vinf_decr":      0.,
    "vrna_decr":      0.,
    "vprot_decr":     0.,
    "die":            0.,
    "k_isg0":         0.,
    "isg_pro_pirf3":  0.,
}


if __name__ == '__main__':
    import click
    import json
    if click.confirm('Save parameters to "./default.json"?'):
        with open('./default.json', 'w') as f:
            json.dump(parameters, f, indent=True)
