# The Dark Matter of Clinical Data — interactive companion

Interactive explainer for the paper *The Dark Matter of Clinical Data: an identification limit on disentangling physiology from care process in electronic health records* (Shahin et al., 2026).

**Live site:** enable GitHub Pages (Settings → Pages → Deploy from branch → `main` / root). The site will appear at `https://<your-username>.github.io/<repo-name>/`.

## What's here
| Path | Contents |
|---|---|
| `index.html` | The interactive page (no build step, no server; all math runs in the browser) |
| `media/` | Narrated explainer videos and posters |
| `paper/` | Main text and Supplementary Information (PDF) |
| `code/` | `reproduce_example.py` and `check_two_tests.py`, which reproduce every number in the paper |

## After the preprint is posted
Open `index.html`, find `const ARXIV_URL = "";` near the top of the script, and paste the arXiv link between the quotes.

## Citation
```bibtex
@misc{shahin2026darkmatter,
  title  = {The Dark Matter of Clinical Data: An identification limit on disentangling physiology from care process in electronic health records},
  author = {Shahin, Mohammad and Chaisutyakorn, Krittaphas and Proios, Dimitrios and Restrepo, David and Mishra, Vishala and Zhang, Boya and Sivaraman, Aadarsh and Roy, Sparsh and Fit{\'e} Abril, Benet and Chen, Tina Yen-Ting and Kwak, Hyunjung Gloria and Kurepa Waschka, Andre and Mathur, Vinayak and Celi, Leo Anthony},
  year   = {2026},
  note   = {Preprint}
}
```
