# Mu2eEvtAna: Analysis using Mu2e Event Ntuple inputs

## Build

```bash
git clone git@github.com:michaelmackenzie/Mu2eEvtAna.git
mu2einit
muse setup
muse build -j<N>
```

## Running

Example processing:
```bash
root.exe -q -b "Mu2eEvtAna/test/test.C"
```