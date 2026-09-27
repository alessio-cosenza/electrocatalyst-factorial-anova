# Publish this project on GitHub from Windows

Suggested repository name: `electrocatalyst-factorial-anova`

Suggested description: `Reproducible Python analysis of catalyst activity using factorial ANOVA, interaction contrasts, and uncertainty estimates.`

## 1. Extract the project

Download `electrocatalyst-factorial-anova.zip` and extract it to your Documents folder. Open PowerShell inside the extracted `electrocatalyst-factorial-anova` directory. Check that `README.md`, `requirements.txt`, `src`, `data`, and `notebooks` are directly inside it.

```powershell
Get-Location
Get-ChildItem
```

If Windows created an extra enclosing folder during extraction, open the inner project folder. Run the analysis using the instructions in the README before publishing.

## 2. Create an empty repository on the intended account

Visit [GitHub's new repository page](https://github.com/new). Check the owner carefully, enter the suggested name and description, and choose **Public**. Leave automatic README, `.gitignore`, and license initialization unchecked because this project already contains files. Click **Create repository** and copy its HTTPS URL.

The correct account is the one whose public profile you want employers to visit. Do not assume a connected application account is the same as your portfolio account.

## 3. Commit the local project

If Git is not installed, install [Git for Windows](https://git-scm.com/download/win), then reopen PowerShell. From the extracted project folder:

```powershell
git init -b main
git add .
git status
git commit -m "Add factorial ANOVA learning project for catalyst activity"
```

If Git asks for your author identity, configure it for this repository using your chosen name and an email associated with your GitHub account. A GitHub-provided no-reply email is suitable if you prefer it. Then repeat the commit command.

```powershell
git config user.name "Alessio Cosenza"
git config user.email "YOUR_GITHUB_EMAIL_OR_NOREPLY_EMAIL"
```

Replace the URL below with the exact URL copied from the new repository:

```powershell
git remote add origin https://github.com/YOUR_USERNAME/electrocatalyst-factorial-anova.git
git push -u origin main
```

Complete GitHub's sign-in flow when prompted. Never paste tokens into tracked files. If the repository was initialized remotely by mistake, stop and reconcile the existing commit history instead of force-pushing.

## Alternative: GitHub CLI

If `gh` is already installed, you can create and push directly after the local commit:

```powershell
gh auth login
gh auth status
gh repo create electrocatalyst-factorial-anova --public --source=. --remote=origin --push --description "Reproducible Python analysis of catalyst activity using factorial ANOVA, interaction contrasts, and uncertainty estimates."
```

Check `gh auth status` before creation to confirm the intended account. This alternative replaces the manual remote-creation steps; do not run both creation routes.

## 4. Check the public result

Open the repository. Confirm that the README image and notebook outputs render, then check the **Actions** tab for the reproducibility workflow. Add relevant topics such as `python`, `anova`, `design-of-experiments`, `statistics`, `materials-science`, and `electrocatalysis`. Pin the repository from your GitHub profile.

Suggested profile description: “Experimental data analysis with Python: factorial modeling, interaction effects, confidence intervals, and transparent scientific reporting.”

## 5. Save future changes

After editing the notebook or analysis code, run all notebook cells, inspect the results, and then:

```powershell
git add .
git diff --cached --stat
git commit -m "Explain factorial effects and update analysis"
git push
```

Reference: [GitHub documentation for adding local code](https://docs.github.com/en/migrations/importing-source-code/using-the-command-line-to-import-source-code/adding-locally-hosted-code-to-github).
