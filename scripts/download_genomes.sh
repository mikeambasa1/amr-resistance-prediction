#!/bin/bash
# Resumable: skips genomes already downloaded; writes to .part so partial files are never mistaken for complete ones.
cd ~/amr-project
mkdir -p data/genomes logs
get() {
  i="$1"
  out="data/genomes/$i.fna"
  [ -s "$out" ] && exit 0
  if curl -s --ssl-reqd --retry 3 --max-time 300 -o "$out.part" "ftp://ftp.bv-brc.org/genomes/$i/$i.fna" && [ -s "$out.part" ]; then
    mv "$out.part" "$out"
  else
    rm -f "$out.part"
    echo "$i" >> logs/failed_downloads.txt
  fi
}
export -f get
cat data/processed/list_complete.txt | xargs -P 4 -I{} bash -c 'get {}'
echo "Done: $(ls data/genomes/*.fna | wc -l) genomes"
