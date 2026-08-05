-- Query for BioSample prefix: SAMEA2
SELECT 
    meta.biosample AS biosample,
    meta.acc AS run_id,
    tax.tax_id AS taxon_id,
    tax.name AS taxon_name,
    tax.self_count AS self_count,
    tax.total_count AS total_count
FROM 
    sra.metadata meta
JOIN 
    sra_tax_analysis_tool.tax_analysis tax ON meta.acc = tax.acc
JOIN
    default.petadex_biosamples pb ON meta.biosample = pb.biosample
WHERE 
    meta.biosample LIKE 'SAMEA2%';
