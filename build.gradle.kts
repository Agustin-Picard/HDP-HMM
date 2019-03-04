import org.jetbrains.kotlin.gradle.tasks.KotlinCompile

plugins {
    java
    kotlin("jvm") version "1.3.21"
}

group = "agustinpicard"
version = "1.0"

repositories {
    mavenCentral()
    maven { url = uri("http://dl.bintray.com/tomasvolker/maven") }
    maven { url = uri("https://jitpack.io") }
    maven {
        url = uri("http://dl.bintray.com/kyonifer/maven")
        jcenter()
    } // koma for extra math functions
}

dependencies {
    implementation(kotlin("stdlib-jdk8"))
    testCompile("junit", "junit", "4.12")

    implementation(group = "tomasvolker", name = "numeriko-core", version = "0.0.3")
    implementation(group = "tomasvolker", name = "kyplot", version = "0.0.1")
    implementation(group = "com.github.tomasvolker", name = "parallel-utils", version = "v1.0")

    compile(group = "com.kyonifer", name = "koma-core-ejml", version = "0.12")
}

configure<JavaPluginConvention> {
    sourceCompatibility = JavaVersion.VERSION_1_8
}
tasks.withType<KotlinCompile> {
    kotlinOptions.jvmTarget = "1.8"
}