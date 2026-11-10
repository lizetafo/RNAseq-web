# Putting the app online (Posit Connect Cloud)

This puts the app on the web so anyone can use it from a link, with no R installation needed. Posit Connect Cloud is free for public GitHub repositories (up to 5 apps, 4 GB memory) and updates the app automatically when you change the code on GitHub.

## 1. Create `manifest.json` (once, on your computer)

Connect Cloud needs a `manifest.json` file that lists the exact R packages the app uses.

1. Open RStudio and set the working directory to this repository folder (**Session → Set Working Directory → Choose Directory…**).
2. Run:
   ```r
   source("setup/make_manifest.R")
   ```
   This installs anything missing, including the six organism databases, which can take 10–20 minutes the first time. It then writes `manifest.json` next to `app.R`.
3. Upload `manifest.json` to the GitHub repository, in the same folder as `app.R`: **Add file → Upload files → Commit changes**.

## 2. Publish

1. Go to **connect.posit.cloud** and sign up for a free account. Signing in with GitHub is easiest.
2. Click **Publish**, then choose **Shiny** as the framework.
3. If asked, install the **Posit Connect Cloud GitHub app** and give it access to `rnaseq-shiny-web`.
4. Choose the repository `lizetafo/rnaseq-shiny-web`, the branch `main`, and `app.R` as the primary file.
5. Click **Publish**. The first build installs all packages and can take 10–30 minutes. Later updates are faster.
6. When it finishes, copy the app's link.

Automatic republishing is on by default, so every time you change the code on GitHub, the online app updates.

## 3. Share the link
- **GitHub:** add the link to the repository's **About → Website** field, and to the top of the README as a **"Try the app online"** link.
- **Website:** add it to the project card on the Projects page as a **"Try the app"** button.
- **LinkedIn:** use it in the Featured section.

## Good to know
- **Data privacy:** files that visitors upload are sent to Posit's servers for that session only. The app does not save them. For unpublished data, run the app on your own computer instead. The app shows this note under the upload buttons.
- **Free-plan resources:** 1 CPU and 4 GB of memory. Datasets with tens of samples work well. Very large studies (hundreds of samples) are better run locally.
- **Sleeping:** free apps may go to sleep when nobody uses them, so the first visitor after a quiet period waits a little while the app wakes up.
- **Try it without data:** the **"Or load example data"** button loads the prostate hypoxia dataset, so visitors can explore the app straight away.

## If something goes wrong
| Problem | What to do |
|---|---|
| `make_manifest.R` stops with a Bioconductor version message | Update BiocManager (`install.packages("BiocManager")`), restart R, run it again. |
| The build fails on Connect Cloud | Open the build **logs** in Connect Cloud. The error near the end usually names the package that failed. Re-run `make_manifest.R` and upload the new `manifest.json`. |
| "Example data not found" | Make sure `example_data/prostate_hypoxia/` is in the repository and listed in `manifest.json` (it is if you used the script). |
