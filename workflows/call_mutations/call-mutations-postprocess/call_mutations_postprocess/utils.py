import collections
import gzip
from pathlib import Path

import bgzip
from click import echo


def fix_clustered_event_flag(
    in_vcf: Path, out_vcf: Path, region_size: int = 100, max_n_somatic_events: int = 2
) -> None:
    """
    Germline variants should not be considered when filtering out clustered somatic
    events. If there are at most `max_n_somatic_events` non-germline events in the
    vicinity of a variant annotated with the `clustered_events` filter (`region_size`),
    then replace the annotation with `PASS`.

    See https://gatk.broadinstitute.org/hc/en-us/community/posts/4404184803227 for more
    information.

    :param in_vcf: a .vcf.gz input file path
    :param out_vcf: a .vcf.gz output file path
    :param max_n_somatic_events: the maximum number of of nearby somatic events to allow
    when potentially repalcing the `clustered_events` flag
    :param region_size: the size of the region in bp to look around
    """

    current_chr = None
    line_buff = collections.deque([])

    with gzip.open(in_vcf, "r+") as f, open(
        out_vcf, "wb"
    ) as raw_out, bgzip.BGZipWriter(raw_out) as f_out:
        for i, line in enumerate(f):
            cur_line = line.decode()

            if cur_line[0] == "#":
                # write header line as-is
                f_out.write(cur_line.encode())
            else:
                # get the current line's columns
                line_cols = cur_line.split("\t")

                if line_cols[0] != current_chr:
                    # starting new chr (or we're starting chr1), so flush the buffer
                    for line_to_write in line_buff:
                        f_out.write("\t".join(line_to_write).encode())

                    # start collecting output lines and keep track of chr and position
                    current_chr = line_cols[0]
                    line_buff = collections.deque([])
                    line_pos = collections.deque([])
                    idx = 0

                    echo(current_chr)

                # continue adding to the buffer
                line_buff.append(line_cols)
                current_pos = int(line_cols[1])
                line_pos.append(current_pos)

                # look ahead until we're more than `region_size` from current position
                while abs(line_pos[idx] - current_pos) > region_size:
                    if line_buff[idx][6] in [
                        "clustered_events",
                        "clustered_events;haplotype",
                    ]:
                        # look around to count other somatic events nearby
                        n_somatic = 0

                        for neighbor in line_buff:
                            if "germline" not in neighbor[6] and abs(
                                int(neighbor[1]) - line_pos[idx]
                            ) <= int(region_size / 2):
                                n_somatic += 1

                        if n_somatic <= max_n_somatic_events:
                            # replace the flag
                            line_buff[idx][6] = "PASS"

                    idx += 1

                # continue to write lines more than 2*region_size from current position
                while abs(line_pos[0] - current_pos) > 2 * region_size:
                    line_to_write = line_buff.popleft()
                    line_cols = "\t".join(line_to_write)
                    f_out.write(line_cols.encode())

                    # roll back the queue
                    idx -= 1
                    line_pos.popleft()

        # write any lines still left in the buffer
        for line_to_write in line_buff:
            f_out.write("\t".join(line_to_write).encode())
