import requests
import json
import time
from dotenv import load_dotenv
import os

def load_api_key():
    load_dotenv()
    api_key = os.getenv("DISGENET_API_KEY")
    if api_key is None:
        raise ValueError("API key not found. Please set the DISGENET_API_KEY environment variable in your .env file.")
    return api_key

def get_gda_summary(api_key, params, max_retries=1):
    headers = {
        "Authorization": api_key,
        "Accept": "application/json"
    }
    url = "https://api.disgenet.com/api/v1/gda/summary"
    for attempt in range(max_retries + 1):
        try:
            response = requests.get(
                url,
                params=params,
                headers=headers,
                timeout=10
            )
        except requests.RequestException as e:
            print(f"Request failed: {e}")
            if attempt < max_retries:
                time.sleep(2)
                continue
            else:
                raise

        if response.status_code == 429 and attempt < max_retries:
            retry_after = int(response.headers.get("x-rate-limit-retry-after-seconds", 60))
            print(f"Rate limit reached, waiting {retry_after}s...")
            time.sleep(retry_after)
            continue
        response.raise_for_status()
        return response.json()
    raise Exception("Failed to get a valid response after retries.")

def save_json(data, filename):
    with open(filename, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    print(f"Summary response saved to {filename}")

def main():
    api_key = load_api_key()
    params = {
        "gene_ncbi_id": "351",
        "page_number": "0"
    }
    response_parsed = get_gda_summary(api_key, params, max_retries=1)
    save_json(response_parsed, "data/summary_response.json")

if __name__ == "__main__":
    main()
