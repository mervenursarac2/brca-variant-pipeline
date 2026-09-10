import requests
from pathlib import Path

URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_PATH = PROJECT_ROOT / "data" / "raw" / "variant_summary.txt.gz"
OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)


def download_with_resume(url, output_path, max_retries=10):
    for attempt in range(max_retries):
        # Ne kadarını zaten indirdiysek, oradan devam et
        existing_size = output_path.stat().st_size if output_path.exists() else 0
        headers = {"Range": f"bytes={existing_size}-"} if existing_size else {}

        try:
            with requests.get(url, headers=headers, stream=True, timeout=30) as r:
                if r.status_code not in (200, 206):
                    print(f"Beklenmeyen durum kodu: {r.status_code}")
                    return False

                total = int(r.headers.get("content-length", 0)) + existing_size
                mode = "ab" if existing_size else "wb"

                with open(output_path, mode) as f:
                    downloaded = existing_size
                    for chunk in r.iter_content(chunk_size=1024 * 1024):  # 1MB parçalar
                        f.write(chunk)
                        downloaded += len(chunk)
                        if total:
                            pct = downloaded / total * 100
                            print(f"\rİndiriliyor: {downloaded/1e6:.1f} MB / {total/1e6:.1f} MB (%{pct:.1f})", end="")

            print("\nİndirme tamamlandı.")
            return True

        except requests.exceptions.RequestException as e:
            print(f"\nDeneme {attempt+1}/{max_retries} kesildi: {e}")
            print("Kaldığı yerden devam edilecek...")

    print("Tüm denemeler tükendi.")
    return False


if __name__ == "__main__":
    download_with_resume(URL, OUTPUT_PATH)