# alpoge-replay-runner

Workflow-only repo for running claimed modules and packets of the volunteer replay
(github.com/rohanarun/alpoge-rh-replay).

1. Open a public claim issue there (`Claim M123456` or `Claim P0123`) and wait for the bot's "Granted" comment.
2. Run `replay-packet` with a single granted ID, or `replay-batch` with several space-separated granted module IDs.
3. Each run uploads its receipts as a workflow artifact, for posting on the claim issue.

No secrets are used and workflows have read-only permissions. `index.html` is published with GitHub Pages
and says who runs this repo. Same-kernel replay is not independent verification.
Not affiliated with the replay's maintainer. No verdict on the underlying claim.
