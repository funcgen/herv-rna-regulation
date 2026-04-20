import csv
import time
from Bio import Entrez
import urllib.error
import http.client
from datetime import datetime

# 1. Authenticate
Entrez.email = "tomas.montserrat@cnag.eu"  # REMEMBER TO CHANGE THIS!

def fetch_abstracts_by_year(output_filename="../results/lncrna_abstracts_dataset_full_no_function_filter.csv", start_year=2005):
    """
    Bypasses the 10k limit by slicing the query by publication year.
    """
    # Broad query: Get ALL literature related to long non-coding RNAs
    query = '("RNA, Long Noncoding"[Mesh] OR "lncRNA"[TIAB] OR "lncRNAs"[TIAB] OR "long noncoding RNA"[TIAB] OR "long non-coding RNA"[TIAB] OR "lincRNA"[TIAB] OR "lincRNAs"[TIAB] OR "ceRNA"[TIAB])'
    current_year = datetime.now().year
    total_downloaded = 0
    
    print(f"Starting year-by-year extraction from {start_year} to {current_year}...\n")

    # Prepare the CSV file
    with open(output_filename, mode='w', newline='', encoding='utf-8') as csv_file:
        fieldnames = ['PMID', 'Year', 'Title', 'Abstract']
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()

        # Loop through each year individually
        for year in range(start_year, current_year + 1):
            print(f"--- Processing Year {year} ---")
            
            # Step A: Search restricted to a specific year using mindate/maxdate
            try:
                search_handle = Entrez.esearch(
                    db="pubmed", 
                    term=query, 
                    usehistory="y",
                    mindate=str(year),
                    maxdate=str(year),
                    datetype="pdat" # Search by Publication Date
                )
                search_results = Entrez.read(search_handle)
                search_handle.close()
            except Exception as e:
                print(f"Error during initial search for {year}: {e}")
                continue

            count = int(search_results["Count"])
            
            if count == 0:
                print(f"No articles found for {year}. Moving on.\n")
                continue
            
            # Safety check just in case a single year explodes in popularity
            if count > 9999:
                print(f"WARNING: {year} has {count} articles. Capping at 9,999 to prevent crash.")
                count = 9999

            webenv = search_results["WebEnv"]
            query_key = search_results["QueryKey"]
            
            print(f"Found {count} articles. Downloading...")
            
            batch_size = 500
            year_downloaded = 0

            # Step B: Batch download for this specific year
            for start in range(0, count, batch_size):
                end = min(count, start + batch_size)
                
                attempt = 0
                success = False
                
                while attempt < 3 and not success:
                    try:
                        fetch_handle = Entrez.efetch(db="pubmed", 
                                                     retmode="xml", 
                                                     retstart=start, 
                                                     retmax=batch_size, 
                                                     webenv=webenv, 
                                                     query_key=query_key)
                        articles_data = Entrez.read(fetch_handle)
                        fetch_handle.close()
                        success = True
                    except Exception as e:
                        attempt += 1
                        print(f"    Network error: HTTP 400 or Timeout. Retrying in 5s (Attempt {attempt}/3)...")
                        time.sleep(5)
                
                if not success:
                    print("    Failed to fetch batch. Skipping to next batch.")
                    continue

                # Step C: Parse and write to CSV
                for article in articles_data.get('PubmedArticle', []):
                    medline = article['MedlineCitation']
                    article_info = medline['Article']
                    
                    pmid = str(medline.get('PMID', 'N/A'))
                    title = article_info.get('ArticleTitle', 'N/A')
                    
                    abstract_text = ""
                    if 'Abstract' in article_info and 'AbstractText' in article_info['Abstract']:
                        abstract_parts = article_info['Abstract']['AbstractText']
                        abstract_text = " ".join([str(part) for part in abstract_parts])
                    
                    if abstract_text:
                        writer.writerow({
                            'PMID': pmid,
                            'Year': str(year),
                            'Title': title,
                            'Abstract': abstract_text
                        })
                        year_downloaded += 1
                
                time.sleep(1) # Polite pause for NCBI

            total_downloaded += year_downloaded
            print(f"Finished {year}. Saved {year_downloaded} valid abstracts.\n")

    print(f"Done! Successfully bypassed the limit. Total abstracts saved: {total_downloaded}")

# --- Run the Script ---
fetch_abstracts_by_year(output_filename="../results/lncrna_abstracts_dataset_full_no_function_filter.csv", start_year=2000)
