vcf <- open_vcf("inst/test.vcf.gz")
vars <- read_vars(vcf, n=10)
vars <- read_vars(vcf, n=10, info_keys=c("AF","AC"), format_keys=c("GT","AD","DP"), format_keys_exclude=c("PL"))
vars <- read_vars(vcf, n=10, info_keys=c("AF","AC"), format_keys=c("GT","AD","DP"), format_keys_exclude=c("PL"), split_info=T)
