from noodles import (gather, schedule)


@schedule
def filter_homo_lumo_lower_than(jobs, x):
    """
    Filter the `jobs` that fulfill that the HOMO-LUMO gap
    is lower than x
    """
    interesting = [(j.job_name, (j.lumo - j.homo) * 27.21) for
                   j in jobs if (j.lumo - j.homo) * 27.21 < x]

    return gather(*interesting)
