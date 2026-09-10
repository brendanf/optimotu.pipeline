test_that("guilds: no and omitted leave guilds off", {
  withr::local_options(list(
    optimotu.pipeline.do_guilds = TRUE,
    optimotu.pipeline.guild_databases = NULL
  ))
  optimotu.pipeline:::parse_guilds_options(list())
  expect_false(optimotu.pipeline::do_guilds())
  expect_equal(nrow(optimotu.pipeline::guild_databases()), 0L)

  optimotu.pipeline:::parse_guilds_options(list(guilds = FALSE))
  expect_false(optimotu.pipeline::do_guilds())
  expect_equal(nrow(optimotu.pipeline::guild_databases()), 0L)
})

test_that("guilds: yes defaults to funguild and optional carlos", {
  withr::local_options(list(
    optimotu.pipeline.do_guilds = FALSE,
    optimotu.pipeline.guild_databases = NULL
  ))
  withr::with_dir(withr::local_tempdir(), {
    expect_message(
      optimotu.pipeline:::parse_guilds_options(list(guilds = TRUE)),
      "skipping that database"
    )
    expect_true(optimotu.pipeline::do_guilds())
    dbs <- optimotu.pipeline::guild_databases()
    expect_equal(dbs$name, "funguild")
    expect_equal(dbs$source, "download")

    dir.create(file.path("data", "lifestyle"), recursive = TRUE)
    saveRDS(
      tibble::tibble(taxon = "Foo", guild = "Pathogen"),
      optimotu.pipeline::lifestyle_guild_path()
    )
    optimotu.pipeline:::parse_guilds_options(list(guilds = TRUE))
    expect_equal(
      optimotu.pipeline::guild_databases()$name,
      c("funguild", "carlos")
    )
    expect_equal(
      optimotu.pipeline::guild_databases()$source,
      c("download", "lifestyle")
    )
  })
})

test_that("guilds list accepts builtins and name/file maps", {
  withr::local_options(list(
    optimotu.pipeline.do_guilds = FALSE,
    optimotu.pipeline.guild_databases = NULL
  ))
  custom <- withr::local_tempfile(fileext = ".tsv")
  readr::write_tsv(
    tibble::tibble(
      taxon = "Foo",
      guild = "Pathogen",
      searchkey = "@Foo@"
    ),
    custom
  )
  withr::with_dir(withr::local_tempdir(), {
    dir.create(file.path("data", "lifestyle"), recursive = TRUE)
    saveRDS(
      tibble::tibble(taxon = "Foo", guild = "Pathogen"),
      optimotu.pipeline::lifestyle_guild_path()
    )
    optimotu.pipeline:::parse_guilds_options(list(
      guilds = list(
        "funguild",
        "carlos",
        list(name = "nemaguild", file = custom)
      )
    ))
  })
  dbs <- optimotu.pipeline::guild_databases()
  expect_equal(dbs$name, c("funguild", "carlos", "nemaguild"))
  expect_equal(dbs$source, c("download", "lifestyle", "file"))
  expect_equal(dbs$path[[3]], custom)
})

test_that("explicit carlos errors when the RDS is missing", {
  withr::with_dir(withr::local_tempdir(), {
    expect_error(
      optimotu.pipeline:::parse_guilds_options(list(guilds = "carlos")),
      "was not found"
    )
  })
})

test_that("custom guild maps reject builtins and missing files", {
  expect_error(
    optimotu.pipeline:::parse_guilds_options(list(
      guilds = list(list(name = "funguild", file = "x.tsv"))
    )),
    "as a string"
  )
  expect_error(
    optimotu.pipeline:::parse_guilds_options(list(
      guilds = list(list(name = "nema", file = "no-such-file.tsv"))
    )),
    "was not found"
  )
  expect_error(
    optimotu.pipeline:::parse_guilds_options(list(guilds = "mfg")),
    "Unknown guild database"
  )
})

test_that("duplicate guild names are rejected", {
  expect_error(
    optimotu.pipeline:::parse_guilds_options(list(
      guilds = c("funguild", "funguild")
    )),
    "Duplicate guild database"
  )
})

test_that("read_guild_table builds searchkey from taxon", {
  path <- withr::local_tempfile(fileext = ".tsv")
  readr::write_tsv(
    tibble::tibble(taxon = "Foo bar", guild = "Pathogen"),
    path
  )
  out <- optimotu.pipeline:::read_guild_table(path)
  expect_equal(out$searchkey, "@Foo@bar@")
})

test_that("guild_db_target_names follows name and file targets", {
  withr::local_options(list(
    optimotu.pipeline.guild_databases = tibble::tibble(
      name = c("funguild", "carlos"),
      source = c("download", "lifestyle"),
      path = c(NA_character_, "data/lifestyle/Fung_LifeStyle_Data.RDS")
    )
  ))
  expect_setequal(
    optimotu.pipeline::guild_db_target_names(),
    c("guild_db_funguild", "guild_db_carlos", "guild_db_file_carlos")
  )
})
