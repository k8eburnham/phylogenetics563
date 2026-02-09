import sys

# The list of 10 species you want to keep
wanted = {
    "Acropyga_sp_FMNHINS0003471605",
    "Anochetus_mayri_FMNHINS0003165063",
    "Camponotus_floridanus_FMNHINS0003145311",
    "Dolichoderus_bispinosus_FMNHINS0003165075",
    "Formica_occulta_FMNHINS0003471213",
    "Iberoformica_subrufa_CASENT0270631",
    "Myrmica_wheeleri_FMNHINS0003471188",
    "Ooceraea_australis_CASENT0106146",
    "Pheidole_flavens_FMNHINS0003145600",
    "Solenopsis_invicta_FMNHINS0003144951"
}

def filter_fasta(filename):
    out_name = "filtered_" + filename
    with open(filename, 'r') as f, open(out_name, 'w') as out:
        keep = False
        for line in f:
            if line.startswith(">"):
                # Extract the name without the '>' and without extra spaces
                header = line[1:].strip().split()[0]
                keep = header in wanted
            if keep:
                out.write(line)
    print(f"Finished {out_name}")

if __name__ == "__main__":
    for arg in sys.argv[1:]:
        filter_fasta(arg)
