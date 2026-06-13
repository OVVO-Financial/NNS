# Downstream sync contract

`OVVO-Financial/NNS` is the source of truth for the NNS implementation. It is
the truth for **both** the native C++ layer under `src/**` and the R API
behavior under `R/**`, `NAMESPACE`, and `DESCRIPTION`.

Downstream flow (a diamond, not a single chain):

```text
OVVO-Financial/NNS  (source of truth: R + C++ src/**)
                   /                       \
        src/** native truth            R API behavior truth
                  v                          v
          NNS-core (portable C++)        (tested directly)
                  \                          |
                   v                         v
                   NNS-python  <-------------+
```

Two paths converge on `NNS-python`:

* **Native path:** `src/**` is the native C++ truth. It flows to `NNS-core`
  (the portable C++ layer), which in turn feeds `NNS-python`.
* **API path:** `R/**`, `NAMESPACE`, and `DESCRIPTION` define the R API
  behavior truth, which is tested directly and flows straight to `NNS-python`
  for API and parity review.

`NNS-python` sits at the convergence of both paths and must satisfy both.

## Triggers

A push to `NNS-Beta-Version` dispatches downstream sync events based on changed files.

| Changed path  | Downstream action                                                             |
| ------------- | ----------------------------------------------------------------------------- |
| `src/**`      | Dispatch `nns-r-src-updated` to `OVVO-Financial/NNS-core`                      |
| `R/**`        | Dispatch `nns-r-api-or-version-updated` to `OVVO-Financial/NNS-python`         |
| `NAMESPACE`   | Dispatch `nns-r-api-or-version-updated` to `OVVO-Financial/NNS-python`         |
| `DESCRIPTION` | Dispatch `nns-r-api-or-version-updated` and require downstream version review  |

A `DESCRIPTION` version change means downstream Python must perform fresh live-R parity-cache regeneration.

## Required secret

`OVVO_SYNC_TOKEN` must be a token or GitHub App installation token with permission to call `repository_dispatch` on:

* `OVVO-Financial/NNS-core`
* `OVVO-Financial/NNS-python`

The default `GITHUB_TOKEN` is intentionally not used for cross-repo dispatch, because it is scoped to this repository and cannot reliably trigger `repository_dispatch` on the downstream repositories.

## Payload

Dispatch payloads include:

* `r_repo`
* `r_commit`
* `r_ref`
* `r_version`
* `r_src_tree_hash`
* `src_changed`
* `r_api_changed`
* `description_changed`
* `changed_files`

## Expected downstream events

* `nns-r-src-updated` — received by `OVVO-Financial/NNS-core` when `src/**` changes.
* `nns-r-api-or-version-updated` — received by `OVVO-Financial/NNS-python` when `R/**`, `NAMESPACE`, or `DESCRIPTION` changes.

## Rule

No downstream sync is complete until each downstream repository records the exact upstream R commit and verifies its own build and parity gates.
