#!/bin/bash
cd ~/amr-project
mkdir -p data/amr logs
JOBS="${1:-2}"
run() {
  i="$1"
  fna="data/genomes/$i.fna"
  out="data/amr/$i.tsv"
  [ -s "$out" ] && exit 0
  [ -s "$fna" ] || exit 0
  if amrfinder -n "$fna" --organism Klebsiella_pneumoniae --plus --threads 1 -o "$out.part" > /dev/null 2>&1 && [ -s "$out.part" ]; then
    mv "$out.part" "$out"
  else
    rm -f "$out.part"
    echo "$i" >> logs/failed_amr.txt
  fi
}
export -f run
cat data/processed/list_complete.txt | xargs -P "$JOBS" -I{} bash -c 'run {}'
echo "Done: $(ls data/amr/*.tsv | wc -l) results"
