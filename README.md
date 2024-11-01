# pyVMO
pyVMO (python Variant Memmap Object), a Python toolkit to help you work with huge variant matrices.


## Installation
```
pip install pyvmo
```

## Usage

### Operation of large vcf files (Recommended use of bcftools to manipulate vcf files)

1. Compress and index your raw vcf files
``` shell
bgzip your_raw.vcf
tabix -p vcf your_raw.vcf.gz
OR
bcftools index --tbi --threads 5 your_raw.vcf.gz
```

2. Extraction of sub-vcf file by a sample list

Only keep the samples in the sample.id.list file and remove sites which are not polymorphic in the sample list.

``` shell
bcftools view --threads 20 -c 1 -O z -o sample.vcf.gz -S sample.id.list raw.vcf.gz
```

3. Extraction of sub-vcf file by variant site quality control

Only keep sites with MAF > 0.05 and missing rate < 0.5.

``` shell
bcftools +fill-tags input.vcf.gz -O z -o input_with_tags.vcf.gz -- -t MAF,F_MISSING
bcftools view -i "MAF>0.05&&F_MISSING<0.5" input_with_tags.vcf.gz -Oz -o filtered_output.vcf.gz
```

4. Extraction of sub-vcf file by region

Only keep sites in the region of chr1:1000000-2000000.

``` shell
bcftools view --threads 20 -c 1 -O z -o region.vcf.gz -r chr1:1000000-2000000 raw.vcf.gz
```

5. Extraction of sub-vcf file by variant site list

Only keep sites in the variant site list file. The variant list file should be a tab-delimited file with the first two columns being the chromosome and position of the variant.

``` shell
bcftools view --threads 20 -c 1 -O z -o varlist.vcf.gz -T var.list raw.vcf.gz
```


### From GWAS

1. Compress and index your raw vcf files
```
bgzip your_raw.vcf
tabix -p vcf your_raw.vcf.gz
```

2. Convert the vcf file to vmo, which is a numpy based memmap file (memory map) that will help you read oversized matrices to get the job done in limited memory

shell
```
PyVMO converter -vcf2vmo your_raw.vcf.gz test.vmo
```

python
```
from pyvmo import VMO

vmo_path = "give_me_your_vmo_path"
raw_vcf_file = "your_raw.vcf.gz"

vmo = VMO(vmo_path)
vmo.store_vcf_to_vmo(raw_vcf_file)
```

3. Extraction of useful submatrices by sample listing, variant site quality control

shell
```
PyVMO extractor sample.id.list test.vmo filter.vmo
```

python
```
# extract by sample list
spl_idx_list = vmo.get_samples_index(sample_id_list)
spl_vmo_path = "give_me_your_sample_list_filtered_vmo_path"
spl_vmo = vmo.extract_sub_vmo(spl_vmo_path, spl_idx_list=spl_idx_list)

# extract by variant site quality control
var_idx_list = spl_vmo.site_filter(maf_thr=0.05, mis_thr=0.5, chunk_size=1000, n_jobs=20)
var_vmo_path = "give_me_your_variant_site_filtered_vmo_path"
var_vmo = spl_vmo.extract_sub_vmo(var_vmo_path, var_idx_list=var_idx_list)
```

4. Convert the vmo file into bimbam format

shell
```
PyVMO converter -vmo2bimbam filter.vmo filter.bimbam
```

python
```
bimbam_file = "give_me_your_bimbam_file_path"
var_vmo.to_bimbam(bimbam_file)
```

5. get genetic distance matrix

shell
```
PyVMO distance filter.vmo ibs.matrix
```

### From Other practices

1. Get numpy array from vmo
```
m = var_vmo.get_matrix()
```

2. Get the sample list
```
sample_list = var_vmo.get_sample_info()
```

3. Get the variant information in a pandas dataframe
```
from yxsql import pickle_load_obj, pickle_dump_obj

var_info_df = var_vmo.get_variant_info()

var_index = 1

chr_id = var_info_df.iloc[i]['CHROM']
pos = int(var_info_df.iloc[i]['POS'])
ref = var_info_df.iloc[i]['REF']
alt = pickle_load_obj(var_info_df.iloc[i]['ALT']) # alt is a list
qual = pickle_load_obj(var_info_df.iloc[i]['QUAL'])
filter = pickle_load_obj(var_info_df.iloc[i]['FILTER'])
info = pickle_load_obj(var_info_df.iloc[i]['INFO']) # info is a dict
```
