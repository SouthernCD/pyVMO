def get_versions():
    return versions[0]["number"]


versions = [
    {
        "number": "0.0.1",
        "features": [
            "1. init",
            "2. VMO builds are included in the initial release, and it is now possible to easily convert vcf files into vmo. with vmo you can easily perform matrix extraction, as well as filter the data according to MAF or Miss ratio.",
            "3. Identity by states (IBS) calculations within a controlled memory and time range",
        ],
    },
]