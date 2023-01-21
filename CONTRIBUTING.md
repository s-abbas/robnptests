# Contributing to `robnptests`
We are grateful for any contribution to `robnptests`.

Before you propose a change via a pull request, you should install the [`devtools`](https://github.com/r-lib/devtools) package as it greatly facilitates the package development process.

## Bugs
If you have found a bug in the code, create a bug report in the [issue tracker](https://github.com/s-abbas/robnptests/issues). Please include a minimal reproducible example ([reprex](https://www.tidyverse.org/help/#reprex)) that illustrates the bug.

## Feature requests
If you miss a feature in the package, just post a new issue in the [issue tracker](https://github.com/s-abbas/robnptests/issues). Before that, please make sure that there is no existing issue with the same request. If such a request already exists, it is, of course, always possible to join the discussion.

Please describe the requested feature as precisely as possible, preferably with a [reprex](https://www.tidyverse.org/help/#reprex).

## Documentation
If you find typos, grammatical errors, spelling or content errors in the documentation, please feel free to make a [pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork).

One important thing to consider is that we use the package [`roxygen2`](https://github.com/r-lib/roxygen2) to document the functions. You should be familiar with the package. 
Changes to the documentation of a function have to be made in the corresponding `.R` file. 
To update the `.Rd` files call the function `devtools::document()`.

## Pull requests
If you would like to actively add a new feature or fix a mistake or bug, you can do so by making a pull request.

For bigger changes, such as extending an existing function or adding a new function, please first create an issue to make sure that someone from the development team agrees that the change is needed.

Follow these steps to contribute to the package:

1. Create a branch in Git from the _develop_ branch with a descriptive name and make your changes.
2. Push the branch to GitHub and create a pull request to merge the changes into the _develop_ branch. Please [link the pull request to the corresponding issue](https://docs.github.com/en/issues/tracking-your-work-with-issues/linking-a-pull-request-to-an-issue).
3. Your pull request will be reviewed by the development team and discussed with you until we approve the pull request or decide that it does not fit for `robnptests`.

Please consider the following points when you create a pull request:

* Pull requests should contain a clear and concise problem description, motivation, and description of the solution.
* Please update the `NEWS.md` file with a description of the change and add your issue number and your GitHub username in the form `(#issue_number, @yourGitHubUsername)` at the end of the description.
* Pull requests should only address one specific change.
  If you would like to contribute multiple changes, please create a new branch for each of it.
  When these changes depend on each other, please submit them separately starting with the first one.
  A new change should only be submitted, after the predecessor has been processed.
* Please make sure that your pull request only contains changes that contribute to the problem the pull request addresses.
* If you change something in the code, it is mandatory that existing unit tests are modified accordingly or new unit tests are added so that the development team and the users of the package can be sure that the new or changed code works as desired.
* New parameters or functions have to be documented with [`roxygen2`](https://github.com/r-lib/roxygen2).
  The function `devtools::document()` has to be called before submitting.
* We strive for a consistent coding style throughout the package.
  Please make sure that your contribution follows the same conventions.
  You do not make anything wrong, if you follow the [tidyverse style guide](https://style.tidyverse.org/).
  If you're modifying existing `robnptests` code that does not follow the style guide, we would appreciate a separate pull request to fix the style.
