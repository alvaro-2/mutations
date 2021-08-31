-- MySQL Script generated by MySQL Workbench
-- Tue Aug 31 16:06:30 2021
-- Model: New Model    Version: 1.0
-- MySQL Workbench Forward Engineering

SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0;
SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0;
SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='ONLY_FULL_GROUP_BY,STRICT_TRANS_TABLES,NO_ZERO_IN_DATE,NO_ZERO_DATE,ERROR_FOR_DIVISION_BY_ZERO,NO_ENGINE_SUBSTITUTION';

-- -----------------------------------------------------
-- Schema mlo_mutations
-- -----------------------------------------------------
DROP SCHEMA IF EXISTS `mlo_mutations` ;

-- -----------------------------------------------------
-- Schema mlo_mutations
-- -----------------------------------------------------
CREATE SCHEMA IF NOT EXISTS `mlo_mutations` DEFAULT CHARACTER SET utf8 ;
USE `mlo_mutations` ;

-- -----------------------------------------------------
-- Table `mlo_mutations`.`protein`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `mlo_mutations`.`protein` ;

CREATE TABLE IF NOT EXISTS `mlo_mutations`.`protein` (
  `id_protein` INT NOT NULL,
  `uniprot_acc` VARCHAR(10) NOT NULL,
  `uniprot_status` VARCHAR(10) NOT NULL,
  `length` INT NULL,
  `uniprot_name` VARCHAR(16) NULL,
  `hgnc_id` VARCHAR(80) NULL,
  `gene_id` VARCHAR(80) NULL,
  `gene_name` VARCHAR(80) NULL,
  `disorder_content` FLOAT NULL,
  `gene_name_synonyms` VARCHAR(300) NULL,
  `protein_names` MEDIUMTEXT NULL,
  `sequence` LONGTEXT NOT NULL,
  PRIMARY KEY (`id_protein`),
  UNIQUE INDEX `uniprot_UNIQUE` (`uniprot_acc` ASC) VISIBLE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `mlo_mutations`.`consequence`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `mlo_mutations`.`consequence` ;

CREATE TABLE IF NOT EXISTS `mlo_mutations`.`consequence` (
  `id_consequence` INT NOT NULL,
  `consequence` VARCHAR(45) NOT NULL,
  PRIMARY KEY (`id_consequence`),
  UNIQUE INDEX `consequence_UNIQUE` (`consequence` ASC) VISIBLE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `mlo_mutations`.`mutation`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `mlo_mutations`.`mutation` ;

CREATE TABLE IF NOT EXISTS `mlo_mutations`.`mutation` (
  `id_mutation` INT NOT NULL,
  `snp_id` VARCHAR(45) NULL,
  `chromosome` VARCHAR(2) NULL,
  `start_genomic` INT NULL,
  `end_genomic` INT NULL,
  `start_aa` INT NULL,
  `end_aa` INT NULL,
  `notation_cds` VARCHAR(300) NULL,
  `notation_aa` VARCHAR(100) NULL,
  `id_protein` INT NOT NULL,
  `id_consequence` INT NOT NULL,
  PRIMARY KEY (`id_mutation`),
  INDEX `fk_mutation_protein_idx` (`id_protein` ASC) VISIBLE,
  INDEX `fk_mutation_consequence_idx` (`id_consequence` ASC) VISIBLE,
  UNIQUE INDEX `notation_UNIQUE` (`notation_aa` ASC, `notation_cds` ASC, `id_protein` ASC) VISIBLE,
  CONSTRAINT `fk_mutation_protein`
    FOREIGN KEY (`id_protein`)
    REFERENCES `mlo_mutations`.`protein` (`id_protein`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_mutation_consequence`
    FOREIGN KEY (`id_consequence`)
    REFERENCES `mlo_mutations`.`consequence` (`id_consequence`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `mlo_mutations`.`disease`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `mlo_mutations`.`disease` ;

CREATE TABLE IF NOT EXISTS `mlo_mutations`.`disease` (
  `id_disease` INT NOT NULL,
  `name` VARCHAR(150) NOT NULL,
  PRIMARY KEY (`id_disease`),
  UNIQUE INDEX `name_UNIQUE` (`name` ASC) VISIBLE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `mlo_mutations`.`disorder_region`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `mlo_mutations`.`disorder_region` ;

CREATE TABLE IF NOT EXISTS `mlo_mutations`.`disorder_region` (
  `id_idr` INT NOT NULL,
  `start` INT NOT NULL,
  `end` INT NOT NULL,
  `length` INT NOT NULL,
  `id_protein` INT NOT NULL,
  PRIMARY KEY (`id_idr`),
  INDEX `fk_disorder_region_protein_idx` (`id_protein` ASC) VISIBLE,
  UNIQUE INDEX `start_end_protein_UNIQUE` (`id_protein` ASC, `start` ASC, `end` ASC) VISIBLE,
  CONSTRAINT `fk_disorder_region_protein`
    FOREIGN KEY (`id_protein`)
    REFERENCES `mlo_mutations`.`protein` (`id_protein`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `mlo_mutations`.`low_complexity`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `mlo_mutations`.`low_complexity` ;

CREATE TABLE IF NOT EXISTS `mlo_mutations`.`low_complexity` (
  `id_lc` INT NOT NULL,
  `start` INT NOT NULL,
  `end` INT NOT NULL,
  `length` INT NOT NULL,
  `id_protein` INT NOT NULL,
  PRIMARY KEY (`id_lc`),
  INDEX `fk_low_complexity_protein_idx` (`id_protein` ASC) INVISIBLE,
  UNIQUE INDEX `start_end_protein_UNIQUE` (`id_protein` ASC, `start` ASC, `end` ASC) VISIBLE,
  CONSTRAINT `fk_low_complexity_protein`
    FOREIGN KEY (`id_protein`)
    REFERENCES `mlo_mutations`.`protein` (`id_protein`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `mlo_mutations`.`pfam_domain`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `mlo_mutations`.`pfam_domain` ;

CREATE TABLE IF NOT EXISTS `mlo_mutations`.`pfam_domain` (
  `id_pfam` VARCHAR(7) NOT NULL,
  `pfam_domain` VARCHAR(20) NOT NULL,
  PRIMARY KEY (`id_pfam`),
  UNIQUE INDEX `pfam_domain_UNIQUE` (`pfam_domain` ASC) VISIBLE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `mlo_mutations`.`mutation_has_disease`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `mlo_mutations`.`mutation_has_disease` ;

CREATE TABLE IF NOT EXISTS `mlo_mutations`.`mutation_has_disease` (
  `id_mutation` INT NOT NULL,
  `id_disease` INT NOT NULL,
  PRIMARY KEY (`id_mutation`, `id_disease`),
  INDEX `fk_mutation_has_disease_disease_idx` (`id_disease` ASC) VISIBLE,
  INDEX `fk_mutation_has_disease_mutation_idx` (`id_mutation` ASC) VISIBLE,
  CONSTRAINT `fk_mutation_has_disease_mutation`
    FOREIGN KEY (`id_mutation`)
    REFERENCES `mlo_mutations`.`mutation` (`id_mutation`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_mutation_has_disease_disease`
    FOREIGN KEY (`id_disease`)
    REFERENCES `mlo_mutations`.`disease` (`id_disease`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `mlo_mutations`.`source`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `mlo_mutations`.`source` ;

CREATE TABLE IF NOT EXISTS `mlo_mutations`.`source` (
  `id_source` INT NOT NULL,
  `source` VARCHAR(45) NOT NULL,
  PRIMARY KEY (`id_source`),
  UNIQUE INDEX `source_UNIQUE` (`source` ASC) VISIBLE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `mlo_mutations`.`mlo`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `mlo_mutations`.`mlo` ;

CREATE TABLE IF NOT EXISTS `mlo_mutations`.`mlo` (
  `id_mlo` INT NOT NULL,
  `mlo` VARCHAR(100) NOT NULL,
  PRIMARY KEY (`id_mlo`),
  UNIQUE INDEX `mlo_UNIQUE` (`mlo` ASC) VISIBLE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `mlo_mutations`.`rol`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `mlo_mutations`.`rol` ;

CREATE TABLE IF NOT EXISTS `mlo_mutations`.`rol` (
  `id_rol` INT NOT NULL,
  `rol` VARCHAR(45) NOT NULL,
  PRIMARY KEY (`id_rol`),
  UNIQUE INDEX `rol_UNIQUE` (`rol` ASC) VISIBLE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `mlo_mutations`.`protein_has_pfam_domain`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `mlo_mutations`.`protein_has_pfam_domain` ;

CREATE TABLE IF NOT EXISTS `mlo_mutations`.`protein_has_pfam_domain` (
  `id_protein` INT NOT NULL,
  `id_pfam` VARCHAR(7) NOT NULL,
  `start` INT NOT NULL,
  `end` INT NOT NULL,
  `length` INT NOT NULL,
  PRIMARY KEY (`id_protein`, `id_pfam`, `start`, `end`),
  INDEX `fk_protein_has_pfam_domain_pfam_domain_idx` (`id_pfam` ASC) VISIBLE,
  INDEX `fk_protein_has_pfam_domain_protein_idx` (`id_protein` ASC) VISIBLE,
  CONSTRAINT `fk_protein_has_pfam_domain_protein`
    FOREIGN KEY (`id_protein`)
    REFERENCES `mlo_mutations`.`protein` (`id_protein`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_protein_has_pfam_domain_pfam_domain`
    FOREIGN KEY (`id_pfam`)
    REFERENCES `mlo_mutations`.`pfam_domain` (`id_pfam`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `mlo_mutations`.`mutation_has_low_complexity`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `mlo_mutations`.`mutation_has_low_complexity` ;

CREATE TABLE IF NOT EXISTS `mlo_mutations`.`mutation_has_low_complexity` (
  `id_mutation` INT NOT NULL,
  `id_lc` INT NOT NULL,
  PRIMARY KEY (`id_mutation`, `id_lc`),
  INDEX `fk_mutation_has_low_complexity_low_complexity_idx` (`id_lc` ASC) VISIBLE,
  INDEX `fk_mutation_has_low_complexity_mutation_idx` (`id_mutation` ASC) VISIBLE,
  CONSTRAINT `fk_mutation_has_low_complexity_mutation`
    FOREIGN KEY (`id_mutation`)
    REFERENCES `mlo_mutations`.`mutation` (`id_mutation`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_mutation_has_low_complexity_low_complexity`
    FOREIGN KEY (`id_lc`)
    REFERENCES `mlo_mutations`.`low_complexity` (`id_lc`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `mlo_mutations`.`mutation_has_disorder_region`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `mlo_mutations`.`mutation_has_disorder_region` ;

CREATE TABLE IF NOT EXISTS `mlo_mutations`.`mutation_has_disorder_region` (
  `id_mutation` INT NOT NULL,
  `id_idr` INT NOT NULL,
  PRIMARY KEY (`id_mutation`, `id_idr`),
  INDEX `fk_mutation_has_disorder_region_disorder_region_idx` (`id_idr` ASC) VISIBLE,
  INDEX `fk_mutation_has_disorder_region_mutation_idx` (`id_mutation` ASC) VISIBLE,
  CONSTRAINT `fk_mutation_has_disorder_region_mutation`
    FOREIGN KEY (`id_mutation`)
    REFERENCES `mlo_mutations`.`mutation` (`id_mutation`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_mutation_has_disorder_region_disorder_region`
    FOREIGN KEY (`id_idr`)
    REFERENCES `mlo_mutations`.`disorder_region` (`id_idr`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `mlo_mutations`.`dataset`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `mlo_mutations`.`dataset` ;

CREATE TABLE IF NOT EXISTS `mlo_mutations`.`dataset` (
  `id_dataset` INT NOT NULL,
  `dataset` VARCHAR(45) NOT NULL,
  PRIMARY KEY (`id_dataset`),
  UNIQUE INDEX `dataset_UNIQUE` (`dataset` ASC) VISIBLE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `mlo_mutations`.`protein_has_mlo`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `mlo_mutations`.`protein_has_mlo` ;

CREATE TABLE IF NOT EXISTS `mlo_mutations`.`protein_has_mlo` (
  `id_proteinmlo` INT NOT NULL,
  `id_protein` INT NOT NULL,
  `id_dataset` INT NOT NULL,
  `id_mlo` INT NULL,
  `id_rol` INT NOT NULL,
  INDEX `fk_protein_has_mlo_mlo_idx` (`id_mlo` ASC) VISIBLE,
  INDEX `fk_protein_has_mlo_protein_idx` (`id_protein` ASC) VISIBLE,
  INDEX `fk_protein_has_mlo_rol_idx` (`id_rol` ASC) VISIBLE,
  INDEX `fk_protein_has_mlo_dataset_idx` (`id_dataset` ASC) VISIBLE,
  UNIQUE INDEX `id_all_UNIQUE` (`id_protein` ASC, `id_dataset` ASC, `id_mlo` ASC) VISIBLE,
  PRIMARY KEY (`id_proteinmlo`),
  CONSTRAINT `fk_protein_has_mlo_protein`
    FOREIGN KEY (`id_protein`)
    REFERENCES `mlo_mutations`.`protein` (`id_protein`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_protein_has_mlo_mlo`
    FOREIGN KEY (`id_mlo`)
    REFERENCES `mlo_mutations`.`mlo` (`id_mlo`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_protein_has_mlo_rol`
    FOREIGN KEY (`id_rol`)
    REFERENCES `mlo_mutations`.`rol` (`id_rol`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_protein_has_mlo_dataset`
    FOREIGN KEY (`id_dataset`)
    REFERENCES `mlo_mutations`.`dataset` (`id_dataset`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `mlo_mutations`.`mutation_has_pfam_domain`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `mlo_mutations`.`mutation_has_pfam_domain` ;

CREATE TABLE IF NOT EXISTS `mlo_mutations`.`mutation_has_pfam_domain` (
  `id_mutation` INT NOT NULL,
  `id_protein` INT NOT NULL,
  `id_pfam` VARCHAR(7) NOT NULL,
  `start` INT NOT NULL,
  `end` INT NOT NULL,
  PRIMARY KEY (`id_mutation`, `id_protein`, `id_pfam`, `start`, `end`),
  INDEX `fk_mutation_has_pfam_domain_protein_has_pfa_idx` (`id_protein` ASC, `id_pfam` ASC, `start` ASC, `end` ASC) VISIBLE,
  INDEX `fk_mutation_has_pfam_domain_mutation_idx` (`id_mutation` ASC) VISIBLE,
  CONSTRAINT `fk_mutation_has_pfam_domain_mutation`
    FOREIGN KEY (`id_mutation`)
    REFERENCES `mlo_mutations`.`mutation` (`id_mutation`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_mutation_has_pfam_domain_protein_has_pfam`
    FOREIGN KEY (`id_protein` , `id_pfam` , `start` , `end`)
    REFERENCES `mlo_mutations`.`protein_has_pfam_domain` (`id_protein` , `id_pfam` , `start` , `end`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `mlo_mutations`.`mutation_has_source`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `mlo_mutations`.`mutation_has_source` ;

CREATE TABLE IF NOT EXISTS `mlo_mutations`.`mutation_has_source` (
  `id_mutation` INT NOT NULL,
  `id_source` INT NOT NULL,
  `id_insource` VARCHAR(15) NOT NULL,
  PRIMARY KEY (`id_mutation`, `id_source`, `id_insource`),
  INDEX `fk_mutation_has_source_source_idx` (`id_source` ASC) VISIBLE,
  INDEX `fk_mutation_has_source_mutation_idx` (`id_mutation` ASC) VISIBLE,
  CONSTRAINT `fk_mutation_has_source_mutation`
    FOREIGN KEY (`id_mutation`)
    REFERENCES `mlo_mutations`.`mutation` (`id_mutation`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_mutation_has_source_source`
    FOREIGN KEY (`id_source`)
    REFERENCES `mlo_mutations`.`source` (`id_source`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `mlo_mutations`.`citation_source`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `mlo_mutations`.`citation_source` ;

CREATE TABLE IF NOT EXISTS `mlo_mutations`.`citation_source` (
  `id_citation_source` INT NOT NULL,
  `name` VARCHAR(45) NULL,
  PRIMARY KEY (`id_citation_source`),
  UNIQUE INDEX `name_UNIQUE` (`name` ASC) VISIBLE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `mlo_mutations`.`mutation_has_citation`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `mlo_mutations`.`mutation_has_citation` ;

CREATE TABLE IF NOT EXISTS `mlo_mutations`.`mutation_has_citation` (
  `id_mutation` INT NOT NULL,
  `id_citation_source` INT NOT NULL,
  `id_citation` VARCHAR(15) NOT NULL,
  PRIMARY KEY (`id_mutation`, `id_citation_source`, `id_citation`),
  INDEX `fk_mutation_has_citation_mutation_idx` (`id_mutation` ASC) VISIBLE,
  INDEX `fk_mutation_has_citation_citation_source_idx` (`id_citation_source` ASC) VISIBLE,
  CONSTRAINT `fk_mutation_has_citation_mutation`
    FOREIGN KEY (`id_mutation`)
    REFERENCES `mlo_mutations`.`mutation` (`id_mutation`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_mutation_has_citation_citation_source`
    FOREIGN KEY (`id_citation_source`)
    REFERENCES `mlo_mutations`.`citation_source` (`id_citation_source`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `mlo_mutations`.`cross_reference`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `mlo_mutations`.`cross_reference` ;

CREATE TABLE IF NOT EXISTS `mlo_mutations`.`cross_reference` (
  `id_cross` INT NOT NULL,
  `name` VARCHAR(45) NOT NULL,
  `url` VARCHAR(50) NULL,
  PRIMARY KEY (`id_cross`),
  UNIQUE INDEX `name_UNIQUE` (`name` ASC) VISIBLE,
  UNIQUE INDEX `url_UNIQUE` (`url` ASC) VISIBLE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `mlo_mutations`.`disease_has_cross_reference`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `mlo_mutations`.`disease_has_cross_reference` ;

CREATE TABLE IF NOT EXISTS `mlo_mutations`.`disease_has_cross_reference` (
  `id_disease` INT NOT NULL,
  `id_cross` INT NOT NULL,
  `id_incross` VARCHAR(20) NOT NULL,
  PRIMARY KEY (`id_disease`, `id_cross`, `id_incross`),
  INDEX `fk_disease_has_cross_reference_cross_reference_idx` (`id_cross` ASC) VISIBLE,
  INDEX `fk_disease_has_cross_reference_disease_idx` (`id_disease` ASC) VISIBLE,
  CONSTRAINT `fk_disease_has_cross_reference_disease`
    FOREIGN KEY (`id_disease`)
    REFERENCES `mlo_mutations`.`disease` (`id_disease`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_disease_has_cross_reference_cross_reference`
    FOREIGN KEY (`id_cross`)
    REFERENCES `mlo_mutations`.`cross_reference` (`id_cross`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `mlo_mutations`.`type_ptm`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `mlo_mutations`.`type_ptm` ;

CREATE TABLE IF NOT EXISTS `mlo_mutations`.`type_ptm` (
  `id_type` INT NOT NULL,
  `type` VARCHAR(20) NULL,
  PRIMARY KEY (`id_type`),
  UNIQUE INDEX `id_type_UNIQUE` (`id_type` ASC) VISIBLE,
  UNIQUE INDEX `type_UNIQUE` (`type` ASC) VISIBLE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `mlo_mutations`.`class_ptm`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `mlo_mutations`.`class_ptm` ;

CREATE TABLE IF NOT EXISTS `mlo_mutations`.`class_ptm` (
  `id_class` INT NOT NULL,
  `class` VARCHAR(20) NULL,
  UNIQUE INDEX `id_class_UNIQUE` (`id_class` ASC) VISIBLE,
  PRIMARY KEY (`id_class`),
  UNIQUE INDEX `class_UNIQUE` (`class` ASC) VISIBLE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `mlo_mutations`.`ptm`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `mlo_mutations`.`ptm` ;

CREATE TABLE IF NOT EXISTS `mlo_mutations`.`ptm` (
  `id_type` INT NOT NULL,
  `id_class` INT NOT NULL,
  `id_protein` INT NOT NULL,
  `pos_aa` INT NOT NULL,
  `mod` VARCHAR(2) NOT NULL,
  `aa` VARCHAR(2) NULL,
  PRIMARY KEY (`id_class`, `id_type`, `id_protein`, `pos_aa`, `mod`),
  INDEX `fk_ptm_type_ptm_idx` (`id_type` ASC) VISIBLE,
  INDEX `fk_ptm_class_ptm_idx` (`id_class` ASC) VISIBLE,
  INDEX `fk_ptm_protein_idx` (`id_protein` ASC) VISIBLE,
  CONSTRAINT `fk_ptm_type_ptm1`
    FOREIGN KEY (`id_type`)
    REFERENCES `mlo_mutations`.`type_ptm` (`id_type`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_ptm_class_ptm1`
    FOREIGN KEY (`id_class`)
    REFERENCES `mlo_mutations`.`class_ptm` (`id_class`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_ptm_protein1`
    FOREIGN KEY (`id_protein`)
    REFERENCES `mlo_mutations`.`protein` (`id_protein`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `mlo_mutations`.`llps_regions`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `mlo_mutations`.`llps_regions` ;

CREATE TABLE IF NOT EXISTS `mlo_mutations`.`llps_regions` (
  `id_llps` INT NOT NULL,
  `start_aa` INT NOT NULL,
  `end_aa` INT NOT NULL,
  `length` INT NULL,
  `id_protein` INT NOT NULL,
  PRIMARY KEY (`id_llps`),
  INDEX `fk_llps_regions_protein_idx` (`id_protein` ASC) INVISIBLE,
  UNIQUE INDEX `start_end_protein_UNIQUE` (`id_protein` ASC, `start_aa` ASC, `end_aa` ASC) VISIBLE,
  CONSTRAINT `fk_llps_regions_protein1`
    FOREIGN KEY (`id_protein`)
    REFERENCES `mlo_mutations`.`protein` (`id_protein`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


SET SQL_MODE=@OLD_SQL_MODE;
SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS;
SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS;
