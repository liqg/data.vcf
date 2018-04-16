library(testthat)

expect(data.vcf:::uniq_char("aadaaffaada","f") == "aadaafaada", "uniq_char for seting char")
expect(data.vcf:::uniq_char("aadaaffaada","") == "adafada", "uniq_char for empty")

