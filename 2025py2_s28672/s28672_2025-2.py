
#!/usr/bin/env python3
from Bio import Entrez, SeqIO
import argparse, os, time, pandas as pd, matplotlib.pyplot as plt

class NCBIRetriever:
    def __init__(self, email, api_key=None):
        Entrez.email = email
        if api_key: Entrez.api_key = api_key
        Entrez.tool = 'BioScript'

    def search_genbank(self, taxid, max_records=100):
        try:
            search_handle = Entrez.esearch(db="nucleotide", term=f"txid{taxid}[Organism]", retmax=max_records)
            search_results = Entrez.read(search_handle)
            search_handle.close()
            return search_results.get("IdList", [])
        except Exception as e:
            print(f"Error searching GenBank: {e}")
            return []

    def fetch_records(self, id_list):
        max_retries = 3
        for attempt in range(max_retries):
            try:
                time.sleep(attempt * 2)
                records = []
                batch_size = min(3, len(id_list))
                for i in range(0, len(id_list), batch_size):
                    batch = id_list[i:i+batch_size]
                    fetch_handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=",".join(batch))
                    batch_records = list(SeqIO.parse(fetch_handle, "genbank"))
                    fetch_handle.close()
                    records.extend(batch_records)
                    time.sleep(1)
                return records
            except Exception as e:
                print(f"Fetch attempt {attempt+1} failed: {e}")
                if attempt == max_retries - 1:
                    return []

    def filter_by_length(self, records, min_length=0, max_length=float('inf')):
        return [r for r in records if min_length <= len(r.seq) <= max_length]

    def generate_csv_report(self, records, output_file):
        try:
            df = pd.DataFrame([{'Accession': r.id, 'Length': len(r.seq), 'Description': r.description} for r in records])
            df.to_csv(output_file, index=False)
            print(f"CSV report saved to {output_file}")
            return df
        except Exception as e:
            print(f"Error generating CSV report: {e}")
            return None

    def visualize_data(self, records, output_file):
        try:
            data = [(r.id, len(r.seq)) for r in records]
            data.sort(key=lambda x: x[1], reverse=True)
            accessions, lengths = zip(*data) if data else ([], [])

            plt.figure(figsize=(10, 5))
            plt.plot(range(len(accessions)), lengths, marker='o')
            plt.xticks(range(len(accessions)), accessions, rotation=90)
            plt.xlabel('GenBank Accession Number')
            plt.ylabel('Sequence Length (bp)')
            plt.title('Sequence Lengths (Longest to Shortest)')
            plt.tight_layout()
            plt.savefig(output_file)
            plt.close()
            print(f"Visualization saved to {output_file}")
            return True
        except Exception as e:
            print(f"Error generating visualization: {e}")
            return False

def main():
    parser = argparse.ArgumentParser(description='Retrieve and analyze genetic sequences from NCBI GenBank.')
    parser.add_argument('--taxid', required=True, help='Taxonomic ID to search for')
    parser.add_argument('--email', required=True, help='Your email for NCBI Entrez')
    parser.add_argument('--api-key', help='NCBI API key (recommended for faster access)')
    parser.add_argument('--min-length', type=int, default=0, help='Minimum sequence length')
    parser.add_argument('--max-length', type=int, default=float('inf'), help='Maximum sequence length')
    parser.add_argument('--max-records', type=int, default=10, help='Maximum number of records to retrieve')
    parser.add_argument('--csv-output', default='genbank_report.csv', help='Output CSV filename')
    parser.add_argument('--chart-output', default='sequence_lengths.png', help='Output chart filename')

    args = parser.parse_args()
    api_key = args.api_key or os.getenv('NCBI_API_KEY')
    retriever = NCBIRetriever(args.email, api_key)

    print(f"Searching GenBank for taxid: {args.taxid}")
    id_list = retriever.search_genbank(args.taxid, args.max_records)
    if not id_list:
        print("No records found.")
        return

    print(f"Found {len(id_list)} records. Fetching data...")
    records = retriever.fetch_records(id_list)
    if not records:
        print("Failed to fetch records.")
        return

    print(f"Retrieved {len(records)} records. Filtering by length...")
    filtered_records = retriever.filter_by_length(records, args.min_length, args.max_length)
    if not filtered_records:
        print(f"No records match the length criteria (min: {args.min_length}, max: {args.max_length}).")
        return

    print(f"{len(filtered_records)} records match the length criteria.")
    retriever.generate_csv_report(filtered_records, args.csv_output)
    retriever.visualize_data(filtered_records, args.chart_output)

if __name__ == "__main__":
    main()
